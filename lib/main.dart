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

/// Clamps a value between [min] and [max].
double clampValue(double x, double min, double max) {
  if (x < min) return min;
  if (x > max) return max;
  return x;
}

/// The main FLIP fluid class (ported from your JavaScript code)
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

// ... (All other simulation methods remain unchanged)
// Ensure that all methods from the user's original FlipFluid class are present here.
}

/// A small holder for scene data
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

  double obstacleX = 1.5; // Initial obstacle position
  double obstacleY = 1.5;
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
  double simHeight = 3.0; // Total domain height as per user's setup
  double cScale = 1.0; // To transform from simulation coords to device coords

  // Track pointer for obstacle
  bool pointerDown = false;

  // Add a flag to prevent multiple initializations
  bool _isSceneInitialized = false;

  @override
  void initState() {
    super.initState();

    // Start the simulation ticker
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
      // Simulation setup is deferred to the build method where canvas size is known
      _isSceneInitialized = true; // Mark as initialized
    }
  }

  /// Initializes the simulation scene based on the canvas size.
  void setupScene(double canvasWidth, double canvasHeight) {
    // Calculate scaling factors
    cScale = canvasHeight / simHeight;
    simWidth = canvasWidth / cScale;

    double tankWidth = simWidth;
    double tankHeight = simHeight;

    double spacing = tankHeight / 100; // Resolution
    double particleRadius = 0.3 * spacing;

    // Initialize FlipFluid
    double density = 1000.0;
    int res = 100;
    double relWaterHeight = 0.4;
    double relWaterWidth = 0.3;

    double h = tankHeight / res;
    double dx = 2.0 * particleRadius;
    double dy = math.sqrt(3.0) / 2.0 * dx;

    int numX = ((relWaterWidth * tankWidth - 2.0 * h - 2.0 * particleRadius) / dx).floor();
    int numY = ((relWaterHeight * tankHeight - 2.0 * h - 2.0 * particleRadius) / dy).floor();
    int maxParticles = numX * numY;

    scene.fluid = FlipFluid(density, tankWidth, tankHeight, h, particleRadius, maxParticles);

    // Fill particle positions
    scene.fluid.numParticles = numX * numY;
    int p = 0;
    for (int i = 0; i < numX; i++) {
      for (int j = 0; j < numY; j++) {
        scene.fluid.particlePos[p++] = h + particleRadius + dx * i + ((j % 2 == 0) ? 0.0 : particleRadius);
        scene.fluid.particlePos[p++] = h + particleRadius + dy * j;
      }
    }

    // Setup grid cells for the tank
    int n = scene.fluid.fNumY;
    for (int i = 0; i < scene.fluid.fNumX; i++) {
      for (int j = 0; j < scene.fluid.fNumY; j++) {
        double val = 1.0; // Fluid by default
        if (i == 0 || i == scene.fluid.fNumX - 1 || j == 0) {
          val = 0.0; // Solid boundary
        }
        scene.fluid.s[i * n + j] = val;
      }
    }

    // Place the obstacle at initial position
    setObstacle(scene.obstacleX, scene.obstacleY, reset: true);
  }

  /// Sets the obstacle's position and updates the simulation grid accordingly.
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
        f.s[i * n + j] = 1.0; // Reset to fluid
        double dx = (i + 0.5) * f.h - x;
        double dy = (j + 0.5) * f.h - y;
        if (dx * dx + dy * dy < r * r) {
          f.s[i * n + j] = 0.0; // Solid obstacle
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

  /// Handles obstacle interaction based on gesture details.
  void _handleObstacleInteraction(DragUpdateDetails details, double cScale, Size canvasSize) {
    // Calculate local position relative to the entire screen
    final RenderBox box = context.findRenderObject() as RenderBox;
    final Offset localPos = box.globalToLocal(details.globalPosition);

    // Adjust for the control panel's width (250 pixels)
    double controlPanelWidth = 250.0;
    double adjustedX = (localPos.dx - controlPanelWidth) / cScale;
    double adjustedY = (canvasSize.height - localPos.dy) / cScale; // Flip Y-axis

    // Clamp to simulation boundaries
    adjustedX = adjustedX.clamp(0.0, simWidth);
    adjustedY = adjustedY.clamp(0.0, simHeight);

    setState(() {
      scene.obstacleX = adjustedX;
      scene.obstacleY = adjustedY;
      // Optionally, set velocities based on drag
    });
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: const Text("FLIP Fluid in Flutter"),
      ),
      body: LayoutBuilder(
        builder: (context, constraints) {
          // Calculate canvas dimensions
          double controlPanelWidth = 250.0;
          double canvasHeight = constraints.maxHeight;
          double canvasWidth = constraints.maxWidth - controlPanelWidth;
          double cScaleLocal = canvasHeight / simHeight;
          double simWidthLocal = canvasWidth / cScaleLocal;

          if (!_isSceneInitialized) {
            setupScene(canvasWidth, canvasHeight);
          }

          return Row(
            children: [
              // 1. Controls Column
              SizedBox(
                width: controlPanelWidth,
                child: SingleChildScrollView(
                  child: Column(
                    children: [
                      // Start/Stop Button
                      ElevatedButton(
                        onPressed: () {
                          setState(() {
                            scene.paused = !scene.paused;
                          });
                        },
                        child: Text(scene.paused ? 'Start' : 'Stop'),
                      ),
                      const SizedBox(height: 10),
                      // Show Particles Switch
                      SwitchListTile(
                        title: const Text("Particles"),
                        value: scene.showParticles,
                        onChanged: (v) {
                          setState(() {
                            scene.showParticles = v;
                          });
                        },
                      ),
                      // Show Grid Switch
                      SwitchListTile(
                        title: const Text("Grid"),
                        value: scene.showGrid,
                        onChanged: (v) {
                          setState(() {
                            scene.showGrid = v;
                          });
                        },
                      ),
                      // Compensate Drift Switch
                      SwitchListTile(
                        title: const Text("Compensate Drift"),
                        value: scene.compensateDrift,
                        onChanged: (v) {
                          setState(() {
                            scene.compensateDrift = v;
                          });
                        },
                      ),
                      // Separate Particles Switch
                      SwitchListTile(
                        title: const Text("Separate Particles"),
                        value: scene.separateParticles,
                        onChanged: (v) {
                          setState(() {
                            scene.separateParticles = v;
                          });
                        },
                      ),
                      // FLIP Ratio Slider
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
              // 2. Canvas Area with GestureDetector
              Expanded(
                child: GestureDetector(
                  onPanDown: (details) {
                    pointerDown = true;
                    _handleObstacleInteraction(details, cScaleLocal, Size(canvasWidth, canvasHeight));
                    setState(() {
                      scene.paused = false; // Unpause when interacting
                    });
                  },
                  onPanUpdate: (details) {
                    if (pointerDown) {
                      _handleObstacleInteraction(details, cScaleLocal, Size(canvasWidth, canvasHeight));
                    }
                  },
                  onPanEnd: (details) {
                    pointerDown = false;
                    scene.obstacleVelX = 0.0;
                    scene.obstacleVelY = 0.0;
                  },
                  child: CustomPaint(
                    painter: FluidPainter(scene, simWidthLocal, simHeight, cScaleLocal),
                    size: Size(canvasWidth, canvasHeight),
                    child: Container(),
                  ),
                ),
              ),
            ],
          );
        },
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
