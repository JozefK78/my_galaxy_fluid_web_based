import 'dart:async';
import 'dart:math' as math;
import 'dart:typed_data';

import 'package:flutter/material.dart';
import 'package:flutter/scheduler.dart'; // For Ticker
import 'dart:ui';



void main() {
  runApp(const MyApp());
}

// --------------------- FLIP FLUID SIMULATION LOGIC ---------------------
enum CellType {
  fluid,
  air,
  solid,
}

double clamp(double x, double min, double max) {
  if (x < min) return min;
  if (x > max) return max;
  return x;
}

// The main FLIP fluid class (ported from your JavaScript code)
class FlipFluid {
  final double density;
  final int fNumX;
  final int fNumY;
  final double h; // Cell size (grid spacing)
  final double fInvSpacing;
  final int fNumCells;

  // Grid-based arrays
  final Float32List u;   // x-velocity
  final Float32List v;   // y-velocity
  final Float32List du;
  final Float32List dv;
  final Float32List prevU;
  final Float32List prevV;
  final Float32List p;   // pressure
  final Float32List s;   // solid mask [1.0 -> fluid, 0.0 -> solid]
  final Int32List cellType;
  final Float32List cellColor;

  // Particle-based arrays
  final int maxParticles;
  final Float32List particlePos;    // x,y
  final Float32List particleColor;  // r,g,b
  final Float32List particleVel;    // vx,vy
  final Float32List particleDensity;
  double particleRestDensity = 0.0;

  // Partitioning for pushing particles apart
  final double particleRadius;
  final double pInvSpacing;
  final int pNumX;
  final int pNumY;
  final int pNumCells;
  final Int32List numCellParticles;
  final Int32List firstCellParticle;
  final Int32List cellParticleIds;

  int numParticles = 0;

  FlipFluid(
      this.density,
      double width,
      double height,
      double spacing,
      this.particleRadius,
      this.maxParticles,
      ) : fNumX = (width / spacing).floor() + 1,
        fNumY = (height / spacing).floor() + 1,
        h = math.max(width / ((width / spacing).floor() + 1),
            height / ((height / spacing).floor() + 1)),
        fInvSpacing = 1.0 / math.max(width / ((width / spacing).floor() + 1),
            height / ((height / spacing).floor() + 1)),
        fNumCells = ((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1),
        u = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        v = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        du = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        dv = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        prevU = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        prevV = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        p = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        s = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        cellType = Int32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        cellColor = Float32List(3 * ((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        particlePos = Float32List(2 * maxParticles),
        particleColor = Float32List(3 * maxParticles),
        particleVel = Float32List(2 * maxParticles),
        particleDensity = Float32List(((width / spacing).floor() + 1) *
            ((height / spacing).floor() + 1)),
        pInvSpacing = 1.0 / (2.2 * particleRadius),
        pNumX = (width * (1.0 / (2.2 * particleRadius))).floor() + 1,
        pNumY = (height * (1.0 / (2.2 * particleRadius))).floor() + 1,
        pNumCells = ((width * (1.0 / (2.2 * particleRadius))).floor() + 1) *
            ((height * (1.0 / (2.2 * particleRadius))).floor() + 1),
        numCellParticles = Int32List(((width * (1.0/(2.2*particleRadius))).floor() + 1) *
            ((height*(1.0/(2.2*particleRadius))).floor() + 1)),
        firstCellParticle = Int32List((((width*(1.0/(2.2*particleRadius))).floor() + 1) *
            ((height*(1.0/(2.2*particleRadius))).floor() + 1)) + 1),
        cellParticleIds = Int32List(maxParticles)
  {
    // Initialize all particle colors to blue
    for (int i = 0; i < maxParticles; i++) {
      particleColor[3 * i + 2] = 1.0;
    }
  }

  void integrateParticles(double dt, double gravity) {
    for (int i = 0; i < numParticles; i++) {
      particleVel[2 * i + 1] += dt * gravity;
      particlePos[2 * i]     += particleVel[2 * i] * dt;
      particlePos[2 * i + 1] += particleVel[2 * i + 1] * dt;
    }
  }

  void pushParticlesApart(int numIters) {
    const double colorDiffusionCoeff = 0.001;
    numCellParticles.fillRange(0, numCellParticles.length, 0);

    // Count how many particles per cell
    for (int i = 0; i < numParticles; i++) {
      double x = particlePos[2 * i];
      double y = particlePos[2 * i + 1];

      int xi = clamp((x * pInvSpacing).floorToDouble(), 0, pNumX - 1).toInt();
      int yi = clamp((y * pInvSpacing).floorToDouble(), 0, pNumY - 1).toInt();
      int cellNr = xi * pNumY + yi;
      numCellParticles[cellNr]++;
    }

    // Prefix sum
    int first = 0;
    for (int i = 0; i < pNumCells; i++) {
      first += numCellParticles[i];
      firstCellParticle[i] = first;
    }
    firstCellParticle[pNumCells] = first; // guard

    // Fill in the cellParticleIds
    for (int i = 0; i < numParticles; i++) {
      double x = particlePos[2 * i];
      double y = particlePos[2 * i + 1];

      int xi = clamp((x * pInvSpacing).floorToDouble(), 0, pNumX - 1).toInt();
      int yi = clamp((y * pInvSpacing).floorToDouble(), 0, pNumY - 1).toInt();
      int cellNr = xi * pNumY + yi;
      firstCellParticle[cellNr]--;
      cellParticleIds[firstCellParticle[cellNr]] = i;
    }

    double minDist = 2.0 * particleRadius;
    double minDist2 = minDist * minDist;

    // The iterative push
    for (int iter = 0; iter < numIters; iter++) {
      for (int i = 0; i < numParticles; i++) {
        double px = particlePos[2 * i];
        double py = particlePos[2 * i + 1];

        int pxi = (px * pInvSpacing).floor();
        int pyi = (py * pInvSpacing).floor();

        int x0 = math.max(pxi - 1, 0);
        int y0 = math.max(pyi - 1, 0);
        int x1 = math.min(pxi + 1, pNumX - 1);
        int y1 = math.min(pyi + 1, pNumY - 1);

        for (int xi = x0; xi <= x1; xi++) {
          for (int yi = y0; yi <= y1; yi++) {
            int cellNr = xi * pNumY + yi;
            int firstIdx = firstCellParticle[cellNr];
            int last = firstCellParticle[cellNr + 1];
            for (int j = firstIdx; j < last; j++) {
              int id = cellParticleIds[j];
              if (id == i) continue;

              double qx = particlePos[2 * id];
              double qy = particlePos[2 * id + 1];
              double dx = qx - px;
              double dy = qy - py;
              double d2 = dx * dx + dy * dy;
              if (d2 > minDist2 || d2 == 0.0) continue;

              double d = math.sqrt(d2);
              double s = 0.5 * (minDist - d) / d;
              dx *= s;
              dy *= s;
              particlePos[2 * i]     -= dx;
              particlePos[2 * i + 1] -= dy;
              particlePos[2 * id]    += dx;
              particlePos[2 * id + 1] += dy;

              // Diffuse colors
              for (int k = 0; k < 3; k++) {
                double color0 = particleColor[3 * i + k];
                double color1 = particleColor[3 * id + k];
                double color = (color0 + color1) * 0.5;
                particleColor[3 * i + k] =
                    color0 + (color - color0) * colorDiffusionCoeff;
                particleColor[3 * id + k] =
                    color1 + (color - color1) * colorDiffusionCoeff;
              }
            }
          }
        }
      }
    }
  }

  void handleParticleCollisions(
      double obstacleX, double obstacleY, double obstacleRadius,
      double obstacleVelX, double obstacleVelY) {
    double h_ = 1.0 / fInvSpacing;
    double r = particleRadius;
    double minDist = obstacleRadius + r;
    double minDist2 = minDist * minDist;

    double minX = h_ + r;
    double maxX = (fNumX - 1) * h_ - r;
    double minY = h_ + r;
    double maxY = (fNumY - 1) * h_ - r;

    for (int i = 0; i < numParticles; i++) {
      double x = particlePos[2 * i];
      double y = particlePos[2 * i + 1];

      double dx = x - obstacleX;
      double dy = y - obstacleY;
      double d2 = dx * dx + dy * dy;

      // Obstacle collision
      if (d2 < minDist2) {
        // Instead of repositioning, we set velocity to that of obstacle
        particleVel[2 * i] = obstacleVelX;
        particleVel[2 * i + 1] = obstacleVelY;
      }

      // Wall collisions
      if (x < minX) {
        x = minX;
        particleVel[2 * i] = 0.0;
      }
      if (x > maxX) {
        x = maxX;
        particleVel[2 * i] = 0.0;
      }
      if (y < minY) {
        y = minY;
        particleVel[2 * i + 1] = 0.0;
      }
      if (y > maxY) {
        y = maxY;
        particleVel[2 * i + 1] = 0.0;
      }
      particlePos[2 * i] = x;
      particlePos[2 * i + 1] = y;
    }
  }

  void updateParticleDensity() {
    int n = fNumY;
    double h_ = h;
    double h1 = fInvSpacing;
    double h2 = 0.5 * h_;

    particleDensity.fillRange(0, particleDensity.length, 0.0);

    for (int i = 0; i < numParticles; i++) {
      double x = particlePos[2 * i];
      double y = particlePos[2 * i + 1];

      x = clamp(x, h_, (fNumX - 1) * h_);
      y = clamp(y, h_, (fNumY - 1) * h_);

      int x0 = ((x - h2) * h1).floor();
      double tx = ((x - h2) - x0 * h_) * h1;
      int x1 = math.min(x0 + 1, fNumX - 2);

      int y0 = ((y - h2) * h1).floor();
      double ty = ((y - h2) - y0 * h_) * h1;
      int y1 = math.min(y0 + 1, fNumY - 2);

      double sx = 1.0 - tx;
      double sy = 1.0 - ty;

      if (x0 >= 0 && y0 >= 0) {
        particleDensity[x0 * n + y0] += sx * sy;
      }
      if (x1 >= 0 && y0 >= 0) {
        particleDensity[x1 * n + y0] += tx * sy;
      }
      if (x1 >= 0 && y1 >= 0) {
        particleDensity[x1 * n + y1] += tx * ty;
      }
      if (x0 >= 0 && y1 >= 0) {
        particleDensity[x0 * n + y1] += sx * ty;
      }
    }

    if (particleRestDensity == 0.0) {
      double sum = 0.0;
      int numFluidCells = 0;
      for (int i = 0; i < fNumCells; i++) {
        if (cellType[i] == CellType.fluid.index) {
          sum += particleDensity[i];
          numFluidCells++;
        }
      }
      if (numFluidCells > 0) {
        particleRestDensity = sum / numFluidCells;
      }
    }
  }

  // toGrid = true => transferring from particles to grid
  // toGrid = false => transferring from grid back to particles
  void transferVelocities(bool toGrid, [double flipRatio = 0.0]) {
    int n = fNumY;
    double h_ = h;
    double h1 = fInvSpacing;
    double h2 = 0.5 * h_;

    if (toGrid) {
      prevU.setAll(0, u);
      prevV.setAll(0, v);

      du.fillRange(0, du.length, 0.0);
      dv.fillRange(0, dv.length, 0.0);
      u.fillRange(0, u.length, 0.0);
      v.fillRange(0, v.length, 0.0);

      // Mark cell types as AIR, then FLUID if a particle is in the cell
      for (int i = 0; i < fNumCells; i++) {
        cellType[i] = (s[i] == 0.0) ? CellType.solid.index : CellType.air.index;
      }
      for (int i = 0; i < numParticles; i++) {
        double x = particlePos[2 * i];
        double y = particlePos[2 * i + 1];

        int xi = clamp((x * h1).floorToDouble(), 0, fNumX - 1).toInt();
        int yi = clamp((y * h1).floorToDouble(), 0, fNumY - 1).toInt();
        int cellNr = xi * n + yi;
        if (cellType[cellNr] == CellType.air.index) {
          cellType[cellNr] = CellType.fluid.index;
        }
      }
    }

    // For each velocity component
    for (int component = 0; component < 2; component++) {
      double dx = (component == 0) ? 0.0 : h2;
      double dy = (component == 0) ? h2 : 0.0;

      Float32List fArr = (component == 0) ? u : v;
      Float32List prevF = (component == 0) ? prevU : prevV;
      Float32List dArr = (component == 0) ? du : dv;

      for (int i = 0; i < numParticles; i++) {
        double x = particlePos[2 * i];
        double y = particlePos[2 * i + 1];

        x = clamp(x, h_, (fNumX - 1) * h_);
        y = clamp(y, h_, (fNumY - 1) * h_);

        int x0 = math.min(((x - dx) * h1).floor(), fNumX - 2);
        double tx = ((x - dx) - x0 * h_) * h1;
        int x1 = math.min(x0 + 1, fNumX - 2);

        int y0 = math.min(((y - dy) * h1).floor(), fNumY - 2);
        double ty = ((y - dy) - y0 * h_) * h1;
        int y1 = math.min(y0 + 1, fNumY - 2);

        double sx = 1.0 - tx;
        double sy = 1.0 - ty;

        double d0 = sx * sy;
        double d1 = tx * sy;
        double d2 = tx * ty;
        double d3 = sx * ty;

        int nr0 = x0 * n + y0;
        int nr1 = x1 * n + y0;
        int nr2 = x1 * n + y1;
        int nr3 = x0 * n + y1;

        if (toGrid) {
          // from particleVel to grid
          double pv = particleVel[2 * i + component];
          fArr[nr0] += pv * d0;
          dArr[nr0] += d0;

          fArr[nr1] += pv * d1;
          dArr[nr1] += d1;

          fArr[nr2] += pv * d2;
          dArr[nr2] += d2;

          fArr[nr3] += pv * d3;
          dArr[nr3] += d3;

        } else {
          // from grid to particleVel (FLIP or PIC)
          int offset = (component == 0) ? n : 1;

          double valid0 = (cellType[nr0] != CellType.air.index ||
              cellType[nr0 - offset] != CellType.air.index)
              ? 1.0
              : 0.0;
          double valid1 = (cellType[nr1] != CellType.air.index ||
              cellType[nr1 - offset] != CellType.air.index)
              ? 1.0
              : 0.0;
          double valid2 = (cellType[nr2] != CellType.air.index ||
              cellType[nr2 - offset] != CellType.air.index)
              ? 1.0
              : 0.0;
          double valid3 = (cellType[nr3] != CellType.air.index ||
              cellType[nr3 - offset] != CellType.air.index)
              ? 1.0
              : 0.0;

          double v_ = particleVel[2 * i + component];
          double d_ = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

          if (d_ > 0.0) {
            double picV = (valid0 * d0 * fArr[nr0] +
                valid1 * d1 * fArr[nr1] +
                valid2 * d2 * fArr[nr2] +
                valid3 * d3 * fArr[nr3]) /
                d_;

            double corr = (valid0 * d0 * (fArr[nr0] - prevF[nr0]) +
                valid1 * d1 * (fArr[nr1] - prevF[nr1]) +
                valid2 * d2 * (fArr[nr2] - prevF[nr2]) +
                valid3 * d3 * (fArr[nr3] - prevF[nr3])) /
                d_;

            double flipV = v_ + corr;
            particleVel[2 * i + component] =
                (1.0 - flipRatio) * picV + flipRatio * flipV;
          }
        }
      }

      if (toGrid) {
        // finalize the grid velocities
        for (int i = 0; i < fArr.length; i++) {
          if (dArr[i] > 0.0) {
            fArr[i] /= dArr[i];
          }
        }

        // restore solid cells
        for (int i = 0; i < fNumX; i++) {
          for (int j = 0; j < fNumY; j++) {
            bool solid = (cellType[i * n + j] == CellType.solid.index);
            if (solid ||
                (i > 0 && cellType[(i - 1) * n + j] == CellType.solid.index)) {
              u[i * n + j] = prevU[i * n + j];
            }
            if (solid ||
                (j > 0 && cellType[i * n + j - 1] == CellType.solid.index)) {
              v[i * n + j] = prevV[i * n + j];
            }
          }
        }
      }
    }
  }

  void solveIncompressibility(
      int numIters, double dt, double overRelaxation, bool compensateDrift) {
    p.fillRange(0, p.length, 0.0);
    prevU.setAll(0, u);
    prevV.setAll(0, v);

    int n = fNumY;
    double cp = density * h / dt;

    for (int iter = 0; iter < numIters; iter++) {
      for (int i = 1; i < fNumX - 1; i++) {
        for (int j = 1; j < fNumY - 1; j++) {
          if (cellType[i * n + j] != CellType.fluid.index) continue;

          int center = i * n + j;
          int left = (i - 1) * n + j;
          int right = (i + 1) * n + j;
          int bottom = i * n + j - 1;
          int top = i * n + j + 1;

          double sx0 = s[left];
          double sx1 = s[right];
          double sy0 = s[bottom];
          double sy1 = s[top];
          double s_ = sx0 + sx1 + sy0 + sy1;
          if (s_ == 0.0) continue;

          double div = u[right] - u[center] + v[top] - v[center];

          if (particleRestDensity > 0.0 && compensateDrift) {
            // Additional correction for compression
            double compression = particleDensity[center] - particleRestDensity;
            if (compression > 0.0) {
              double k = 1.0;
              div = div - k * compression;
            }
          }

          double p_ = -div / s_;
          p_ *= overRelaxation;
          p[center] += cp * p_;

          u[center] -= sx0 * p_;
          u[right]  += sx1 * p_;
          v[center] -= sy0 * p_;
          v[top]    += sy1 * p_;
        }
      }
    }
  }

  void updateParticleColors() {
    double h1 = fInvSpacing;

    for (int i = 0; i < numParticles; i++) {
      // fade to white (example logic)
      double s = 0.01;
      particleColor[3 * i]     = clamp(particleColor[3 * i] - s, 0.0, 1.0);
      particleColor[3 * i + 1] = clamp(particleColor[3 * i + 1] - s, 0.0, 1.0);
      particleColor[3 * i + 2] = clamp(particleColor[3 * i + 2] + s, 0.0, 1.0);

      double x = particlePos[2 * i];
      double y = particlePos[2 * i + 1];
      int xi = clamp((x * h1).floorToDouble(), 1, fNumX - 1).toInt();
      int yi = clamp((y * h1).floorToDouble(), 1, fNumY - 1).toInt();
      int cellNr = xi * fNumY + yi;

      double d0 = particleRestDensity;
      if (d0 > 0.0) {
        double relDensity = particleDensity[cellNr] / d0;
        if (relDensity < 0.7) {
          double s_ = 0.8;
          particleColor[3 * i]     = s_;
          particleColor[3 * i + 1] = s_;
          particleColor[3 * i + 2] = 1.0;
        }
      }
    }
  }

  void setSciColor(int cellNr, double val, double minVal, double maxVal) {
    val = math.min(math.max(val, minVal), maxVal - 0.0001);
    double d = maxVal - minVal;
    double norm = (d == 0.0) ? 0.5 : (val - minVal) / d;
    double m = 0.25;
    int num = (norm / m).floor();
    double s = (norm - num * m) / m;
    double r, g, b;
    switch (num) {
      case 0:
        r = 0.0;
        g = s;
        b = 1.0;
        break;
      case 1:
        r = 0.0;
        g = 1.0;
        b = 1.0 - s;
        break;
      case 2:
        r = s;
        g = 1.0;
        b = 0.0;
        break;
      case 3:
      default:
        r = 1.0;
        g = 1.0 - s;
        b = 0.0;
        break;
    }
    cellColor[3 * cellNr]     = r;
    cellColor[3 * cellNr + 1] = g;
    cellColor[3 * cellNr + 2] = b;
  }

  void updateCellColors() {
    // Clear
    cellColor.fillRange(0, cellColor.length, 0.0);
    for (int i = 0; i < fNumCells; i++) {
      if (cellType[i] == CellType.solid.index) {
        cellColor[3 * i]     = 0.5;
        cellColor[3 * i + 1] = 0.5;
        cellColor[3 * i + 2] = 0.5;
      } else if (cellType[i] == CellType.fluid.index) {
        double d = particleDensity[i];
        if (particleRestDensity > 0.0) d /= particleRestDensity;
        setSciColor(i, d, 0.0, 2.0);
      }
    }
  }

  // The main time-step
  void simulate(
      double dt,
      double gravity,
      double flipRatio,
      int numPressureIters,
      int numParticleIters,
      double overRelaxation,
      bool compensateDrift,
      bool separateParticles,
      double obstacleX,
      double obstacleY,
      double obstacleRadius,
      double obstacleVelX,
      double obstacleVelY) {
    int numSubSteps = 1;
    double sdt = dt / numSubSteps;
    for (int step = 0; step < numSubSteps; step++) {
      integrateParticles(sdt, gravity);
      if (separateParticles) {
        pushParticlesApart(numParticleIters);
      }
      handleParticleCollisions(obstacleX, obstacleY, obstacleRadius,
          obstacleVelX, obstacleVelY);
      transferVelocities(true);
      updateParticleDensity();
      solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
      transferVelocities(false, flipRatio);
    }
    updateParticleColors();
    updateCellColors();
  }
}

// A small holder for scene data
class Scene {
  double gravity = -9.81;
  double dt = 1.0 / 60.0;
  double flipRatio = 0.9;
  int numPressureIters = 50;
  int numParticleIters = 2;
  double overRelaxation = 1.9;
  bool compensateDrift = true;
  bool separateParticles = true;
  bool paused = true;

  bool showParticles = true;
  bool showGrid = false;

  double obstacleX = 0.0;
  double obstacleY = 0.0;
  double obstacleRadius = 0.15;
  double obstacleVelX = 0.0;
  double obstacleVelY = 0.0;

  // The fluid
  late FlipFluid fluid;
}

// --------------------- FLUTTER APP UI & RENDERING ---------------------

class MyApp extends StatelessWidget {
  const MyApp({super.key});

  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      title: 'FLIP Fluid Demo',
      theme: ThemeData(
        primarySwatch: Colors.blue,
      ),
      home: const FluidDemoPage(),
    );
  }
}

class FluidDemoPage extends StatefulWidget {
  const FluidDemoPage({Key? key}) : super(key: key);

  @override
  State<FluidDemoPage> createState() => _FluidDemoPageState();
}

class _FluidDemoPageState extends State<FluidDemoPage>
    with SingleTickerProviderStateMixin {
  final Scene scene = Scene();
  late Ticker _ticker;

  // Canvas scale factors
  double simWidth = 0;
  double simHeight = 0;
  double cScale = 1.0; // to transform from sim coords to device coords

  // Track pointer for obstacle
  bool pointerDown = false;
  double lastPointerX = 0.0;
  double lastPointerY = 0.0;

  // Add a flag to prevent multiple initializations
  bool _isSceneInitialized = false;

  @override
  void initState() {
    super.initState();

    // Remove setupScene() from initState()
    // setupScene();

    // Create a Ticker to update at ~60 fps
    _ticker = createTicker((elapsed) {
      if (!scene.paused) {
        scene.fluid.simulate(
          scene.dt,
          scene.gravity,
          scene.flipRatio,
          scene.numPressureIters,
          scene.numParticleIters,
          scene.overRelaxation,
          scene.compensateDrift,
          scene.separateParticles,
          scene.obstacleX,
          scene.obstacleY,
          scene.obstacleRadius,
          scene.obstacleVelX,
          scene.obstacleVelY,
        );
        setState(() {});
      }
    });
    _ticker.start();
  }

  @override
  void dispose() {
    _ticker.dispose();
    super.dispose();
  }

  @override
  void didChangeDependencies() {
    super.didChangeDependencies();

    if (!_isSceneInitialized) {
      simHeight = 3.0; // total domain height
      simWidth = (MediaQuery.of(context).size.width /
          MediaQuery.of(context).size.height) *
          simHeight; // keep aspect ratio
      setupScene();
      _isSceneInitialized = true; // Mark as initialized
    }
  }

  void setupScene() {
    // Similar to the JS code:

    scene.obstacleRadius = 0.15;
    scene.overRelaxation = 1.9;
    scene.dt = 1.0 / 60.0;
    scene.numPressureIters = 50;
    scene.numParticleIters = 2;

    int res = 100; // number of cells in the vertical direction
    double tankHeight = simHeight;
    double tankWidth = simWidth;

    double h = tankHeight / res;
    double density = 1000.0;

    double relWaterHeight = 0.4;
    double relWaterWidth = 0.3;

    double r = 0.3 * h; // particle radius
    double dx = 2.0 * r;
    double dy = math.sqrt(3.0) / 2.0 * dx;

    int numX = ((relWaterWidth * tankWidth - 2.0 * h - 2.0 * r) / dx).floor();
    int numY = ((relWaterHeight * tankHeight - 2.0 * h - 2.0 * r) / dy).floor();
    int maxParticles = numX * numY;

    scene.fluid = FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

    // Fill particle positions
    scene.fluid.numParticles = numX * numY;
    int p = 0;
    for (int i = 0; i < numX; i++) {
      for (int j = 0; j < numY; j++) {
        scene.fluid.particlePos[p++] = h + r + dx * i + ((j % 2 == 0) ? 0.0 : r);
        scene.fluid.particlePos[p++] = h + r + dy * j;
      }
    }

    // Setup grid cells for the tank
    int n = scene.fluid.fNumY;
    for (int i = 0; i < scene.fluid.fNumX; i++) {
      for (int j = 0; j < scene.fluid.fNumY; j++) {
        double val = 1.0; // fluid by default
        if (i == 0 || i == scene.fluid.fNumX - 1 || j == 0) {
          val = 0.0; // solid boundary
        }
        scene.fluid.s[i * n + j] = val;
      }
    }

    // Place the obstacle somewhere
    setObstacle(scene.obstacleX, scene.obstacleY, reset: true);

    // We'll compute a scale so that we can fill the screen
    double screenWidth = MediaQuery.of(context).size.width;
    double screenHeight = MediaQuery.of(context).size.height;
    cScale = screenHeight / simHeight; // maintain aspect ratio
  }

  void setObstacle(double x, double y, {bool reset = false}) {
    double vx = 0.0;
    double vy = 0.0;
    if (!reset) {
      vx = (x - scene.obstacleX) / scene.dt;
      vy = (y - scene.obstacleY) / scene.dt;
    }
    scene.obstacleX = x;
    scene.obstacleY = y;
    double r = scene.obstacleRadius;
    var f = scene.fluid;
    int n = f.fNumY;

    for (int i = 1; i < f.fNumX - 2; i++) {
      for (int j = 1; j < f.fNumY - 2; j++) {
        f.s[i * n + j] = 1.0;
        double dx = (i + 0.5) * f.h - x;
        double dy = (j + 0.5) * f.h - y;
        if (dx * dx + dy * dy < r * r) {
          f.s[i * n + j] = 0.0;
          f.u[i * n + j] = vx;
          f.u[(i + 1) * n + j] = vx;
          f.v[i * n + j] = vy;
          f.v[i * n + j + 1] = vy;
        }
      }
    }
    scene.obstacleVelX = vx;
    scene.obstacleVelY = vy;
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text("FLIP Fluid in Flutter"),
      ),
      body: GestureDetector(
        onPanDown: (details) {
          pointerDown = true;
          final localPos = details.localPosition;
          lastPointerX = localPos.dx;
          lastPointerY = localPos.dy;
          double simX = localPos.dx / cScale;
          double simY = (MediaQuery.of(context).size.height - localPos.dy) / cScale;
          setObstacle(simX, simY, reset: true);
          setState(() {
            scene.paused = false; // unpause
          });
        },
        onPanUpdate: (details) {
          if (pointerDown) {
            final localPos = details.localPosition;
            double simX = localPos.dx / cScale;
            double simY = (MediaQuery.of(context).size.height - localPos.dy) / cScale;
            setObstacle(simX, simY, reset: false);
            setState(() {});
          }
        },
        onPanEnd: (details) {
          pointerDown = false;
          scene.obstacleVelX = 0.0;
          scene.obstacleVelY = 0.0;
        },
        child: Row(
          children: [
            // Left column with controls
            SizedBox(
              width: 250,
              child: SingleChildScrollView(
                child: Column(
                  children: [
                    // Start/Stop
                    ElevatedButton(
                      onPressed: () {
                        setState(() {
                          scene.paused = !scene.paused;
                        });
                      },
                      child: Text(scene.paused ? 'Start' : 'Stop'),
                    ),
                    const SizedBox(height: 10),
                    // Show Particles
                    SwitchListTile(
                      title: const Text("Particles"),
                      value: scene.showParticles,
                      onChanged: (v) {
                        setState(() {
                          scene.showParticles = v;
                        });
                      },
                    ),
                    // Show Grid
                    SwitchListTile(
                      title: const Text("Grid"),
                      value: scene.showGrid,
                      onChanged: (v) {
                        setState(() {
                          scene.showGrid = v;
                        });
                      },
                    ),
                    // Compensate Drift
                    SwitchListTile(
                      title: const Text("Compensate Drift"),
                      value: scene.compensateDrift,
                      onChanged: (v) {
                        setState(() {
                          scene.compensateDrift = v;
                        });
                      },
                    ),
                    // Separate Particles
                    SwitchListTile(
                      title: const Text("Separate Particles"),
                      value: scene.separateParticles,
                      onChanged: (v) {
                        setState(() {
                          scene.separateParticles = v;
                        });
                      },
                    ),
                    // PIC -> FLIP ratio
                    const Text("FLIP Ratio"),
                    Slider(
                      min: 0.0,
                      max: 1.0,
                      divisions: 10,
                      value: scene.flipRatio,
                      onChanged: (val) {
                        setState(() {
                          scene.flipRatio = val;
                        });
                      },
                    ),
                  ],
                ),
              ),
            ),
            // The canvas (CustomPaint)
            Expanded(
              child: CustomPaint(
                painter: FluidPainter(scene, simWidth, simHeight, cScale),
                child: Container(),
              ),
            ),
          ],
        ),
      ),
    );
  }
}

// -------------------- Rendering via CustomPainter --------------------
class FluidPainter extends CustomPainter {
  final Scene scene;
  final double simWidth;
  final double simHeight;
  final double cScale;

  // Reusable Paint objects to optimize performance
  final Paint _backgroundPaint = Paint()..color = Colors.black;
  final Paint _gridPaint = Paint();
  final Paint _particlePaint = Paint()
    ..color = Colors.white
    ..style = PaintingStyle.fill
    ..strokeCap = StrokeCap.round
    ..strokeWidth = 2.0; // Adjust based on particle size

  final Paint _obstaclePaint = Paint()
    ..color = Colors.red
    ..style = PaintingStyle.fill;

  FluidPainter(this.scene, this.simWidth, this.simHeight, this.cScale);

  @override
  void paint(Canvas canvas, Size size) {
    // 1. Draw Background
    canvas.drawRect(
      Rect.fromLTWH(0, 0, size.width, size.height),
      _backgroundPaint,
    );

    final FlipFluid fluid = scene.fluid;

    // 2. Draw Grid (Optional)
    if (scene.showGrid) {
      _drawGrid(canvas, fluid);
    }

    // 3. Draw Particles (Optional)
    if (scene.showParticles) {
      _drawParticles(canvas, fluid);
    }

    // 4. Draw Obstacle
    _drawObstacle(canvas, fluid);
  }

  /// Draws the simulation grid.
  void _drawGrid(Canvas canvas, FlipFluid f) {
    final int n = f.fNumY; // Number of cells vertically

    for (int i = 0; i < f.fNumX; i++) {
      for (int j = 0; j < f.fNumY; j++) {
        int cellNr = i * n + j;

        // Retrieve color components for the cell
        double r = f.cellColor[3 * cellNr];
        double g = f.cellColor[3 * cellNr + 1];
        double b = f.cellColor[3 * cellNr + 2];

        // Set the paint color
        _gridPaint.color = Color.fromARGB(
          255,
          (r * 255).floor(),
          (g * 255).floor(),
          (b * 255).floor(),
        );

        // Calculate the position of the cell
        double x = i * f.h;
        double y = j * f.h;
        double sy = simHeight - (y + f.h); // Flip Y-axis for Flutter

        // Define the rectangle for the cell
        final Rect rect = Rect.fromLTWH(
          x * cScale,
          sy * cScale,
          f.h * cScale,
          f.h * cScale,
        );

        // Draw the cell
        canvas.drawRect(rect, _gridPaint);
      }
    }
  }

  /// Draws all the particles using a single drawPoints call for efficiency.
  void _drawParticles(Canvas canvas, FlipFluid f) {
    final List<Offset> points = <Offset>[];

    for (int i = 0; i < f.numParticles; i++) {
      double x = f.particlePos[2 * i];
      double y = f.particlePos[2 * i + 1];
      double sy = simHeight - y; // Flip Y-axis

      // Ensure particles are within simulation bounds
      if (x < 0 || x > simWidth || y < 0 || y > simHeight) continue;

      points.add(Offset(x * cScale, sy * cScale));
    }

    // Draw all particles in a single draw call
    canvas.drawPoints(PointMode.points, points, _particlePaint);
  }

  /// Draws the obstacle within the simulation.
  void _drawObstacle(Canvas canvas, FlipFluid f) {
    double x = scene.obstacleX * cScale;
    double y = (simHeight - scene.obstacleY) * cScale; // Flip Y-axis
    double radius = (scene.obstacleRadius + f.particleRadius) * cScale;

    // Draw the obstacle as a circle
    canvas.drawCircle(
      Offset(x, y),
      radius,
      _obstaclePaint,
    );
  }

  @override
  bool shouldRepaint(covariant FluidPainter oldDelegate) {
    return true; // Repaint every frame for continuous animation
  }
}
