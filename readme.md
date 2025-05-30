# Pinto: Tree-Jumping Squirrel Robot

![Pinto Robot Image Placeholder](readme_media/treejumpclipi2.gif)

## Overview

Pinto is a 450g tree-jumping robot inspired by squirrel locomotion, developed to demonstrate the feasibility of arboreal robotics. The robot features a novel latched series-elastic actuator mechanism that can switch between stiff series elastic and parallel elastic modes, enabling both precision movement and powerful energy storage for jumping from the ground onto vertical tree trunks.

**Key Features:**
- **Weight:** 450g lightweight design
- **Novel Mechanism:** Latched series-elastic actuator using twisted string and carbon fiber springs
- **5-bar leg mechanism:** Front and rear legs with 2 motors cooperatively actuating each foot
- **Dual-mode operation:** Series and parallel-elastic modes for different locomotion requirements
- **Spined grippers:** 2-DoF arms for grasping tree bark during high-speed perching

**Research Abstract:** Arboreal environments challenge current robots but are deftly traversed by many familiar animal locomotors such as squirrels. We present a small, 450 g robot "Pinto" developed for tree-jumping, a behavior seen in squirrels but rarely in legged robots: jumping from the ground onto a vertical tree trunk. We develop a powerful and lightweight latched series-elastic actuator using a twisted string and carbon fiber springs. We consider the effects of scaling down conventional quadrupeds and experimentally show how storing energy in a parallel-elastic fashion using a latch increases jump energy compared to series-elastic or springless strategies. By switching between series and parallel-elastic modes with our latched 5-bar leg mechanism, Pinto executes energetic jumps as well as maintains continuous control during shorter bounding motions. We also develop sprung 2-DoF arms equipped with spined grippers to grasp tree bark for high-speed perching following a jump.

## Repository Structure

### `pintocode/`
Contains firmware and control code for various microcontroller platforms and iterations:

#### **Version 1: Teensy 4.0 with Thumbstick Control**
- **`teensy/thumbstick*/`**: Early single-leg prototype controlled by thumbstick
- **Platform:** Teensy 4.0 microcontroller
- **Files:** `.ino` files for Arduino IDE

![Version 1 Image Placeholder](placeholder_v1_image.jpg)
*[Placeholder for Version 1 robot image]*

#### **Version 2: ESP8266 Wireless Control System**
- **`teensy/joystick*/`**: Teensy 4.0 code for custom 8-axis joystick (four 2-axis thumbsticks)
- **`esp8266/espnow_send/`**: ESP8266 transceiver code for control station
- **`esp8266/espnow_recv/`**: ESP8266 receiver code for robot
- **Communication:** ESP-NOW wireless protocol
- **Control:** Human input via 8-axis joystick mapped to 4 robot limbs

![Version 2 Image Placeholder](placeholder_v2_image.jpg)
*[Placeholder for Version 2 system image]*

#### **Version 3: Station-Estop-Brain Architecture (Current)**
- **`station/`**: Python control software running on laptop
  - `station.py`: Main control program with gamepad input via pygame
  - `periodics.py`: Periodic task scheduling utilities
  - `data/`: Logged experimental data and trajectories
- **`esp32s3/estop/`**: Emergency stop device on Seeedstudio XIAO ESP32S3
- **`esp32s3/brain2/`**: Main robot brain on custom "squirrelbrain" PCB
- **Communication Chain:** Laptop ↔ USB Serial ↔ Estop ↔ ESP-NOW ↔ Brain

**System Components:**
1. **Station (Laptop):** Reads Bluetooth gamepad, computes joint trajectories and motor commands
2. **Estop (ESP32S3):** Physical emergency stop button with wireless transceiver functionality
3. **Brain (ESP32S3):** Custom PCB controlling servos, brushless motors, and sensors

![Version 3 Station Image Placeholder](placeholder_v3_station.jpg)
*[Placeholder for Version 3 station setup image]*

![Version 3 Estop Image Placeholder](placeholder_v3_estop.jpg)
*[Placeholder for Version 3 estop device image]*

![Version 3 Brain Image Placeholder](placeholder_v3_brain.jpg)
*[Placeholder for Version 3 brain PCB image]*

#### **Other Platforms:**
- **`esp32s3/`**: Various ESP32S3 applications and utilities
- **`seeeduino/`**: Seeeduino-specific firmware
- **`lookuptables/`**: Precomputed lookup tables for inverse kinematics and control

### `pintomech/`
Contains mechanical analysis, simulation, and experimental data:

#### **Mechanism Analysis:**
- **`SEA/`**: Series Elastic Actuator analysis and optimization notebooks
- **`spring_theory/`**: Theoretical spring mechanism analysis
- **`spring_experiments/`**: Physical spring testing and validation data
- **`string_experiments/`**: Twisted string actuator characterization
- **`lambda_linkage/`**: Lambda-linkage mechanism analysis for a vertically moving boom for testing a jumping leg
- **`IK/`**: Inverse kinematics implementations and verification

#### **Research Data:**
- **`ICRA24/`**: 2024 Conference data - force sensor jump experiments
- **`ICRA25/`**: 2025 Conference data and paper materials
  - `spring_torque/`: Torsion spring torque-displacement profiles
  - `cycle/`: Rear foot cyclic tracking experiments
  - `string_mapping/`: Twisted string input-output characteristics
  - `trials_used/`: Jump experiments with camera tracking
- **`jump_force/`**: Jump force analysis and experiments
- **`polytope/`**: Polytope-based motion planning analysis
- **`windup/`**: Windup mechanism studies

![Jump Experiment Image Placeholder](placeholder_jump_experiment.jpg)
*[Placeholder for jump experiment setup image]*

#### **Utilities:**
- **`logs/`**: Experimental logs and data files

## Hardware Requirements

- **Custom PCBs:** Available at [github.com/qwertpas/squirrelbrain](https://github.com/qwertpas/squirrelbrain)
- **CAD Files:** Available at [Onshape CAD](https://cad.onshape.com/documents/0c72dae6c9475dd41cab7700/w/2583ef1f629f03790729c107/e/f951a8c75300fed560aa6177)
- **Microcontrollers:** Teensy 4.0, ESP8266, ESP32S3 (Seeedstudio XIAO)

### **Motors:**
- **Brushless Motors:** [Placeholder for BLDC motor specifications]
- **Servo Motors:** [Placeholder for servo motor specifications - Dynamixel compatible]

### **Sensors:**
- **Encoders:** [Placeholder for encoder specifications]
- **IMU:** [Placeholder for IMU specifications]
- **Force Sensors:** [Placeholder for force sensor specifications]

### **Power:**
- **Battery:** [Placeholder for battery specifications - 3S LiPo, >11V required]
- **Power Management:** [Placeholder for power distribution details]

### **Additional Components:**
- **Springs:** Custom carbon fiber springs
- **Actuators:** Twisted string actuators
- **Grippers:** Spined grippers for bark grasping
- **Controller:** BINBOK Nintendo Switch compatible controller (or Xbox controller)

## Software Requirements

### Python Dependencies (Station)
**Python Version:** 3.11

Install required packages:
```bash
pip install pygame pandas numpy pyserial
```

### Microcontroller Development
- **PlatformIO:** For ESP32S3 and ESP8266 development
- **Arduino IDE:** For legacy Teensy `.ino` files
- **Setup Guide:** [PlatformIO VSCode Integration](https://docs.platformio.org/en/latest/integration/ide/vscode.html)

Dependencies are specified in the `lib_deps` section of each project's `platformio.ini` file.

## Operating Procedure (Latest version)

**Required Code:**
- **Brain:** `pintocode/esp32s3/brain2/src/main.cpp`
- **Station:** `pintocode/station/station.py`
- **Estop:** `pintocode/esp32s3/estop/src/main.cpp`

**Setup Steps:**

1. **Connect Estop to Laptop**
   - Connect estop device to laptop via USB
   - Determine COM port (Windows) or serial port (Linux/Mac)
   - Update the `port` variable in `station.py`
   - Verify small yellow light near USB port is blinking on estop

2. **Run Station Software**
   - Execute `python station.py`
   - Should display: "estop_ok" and "gamepad disconnected"

3. **Connect Gamepad**
   - Connect BINBOK Nintendo Switch controller (or Xbox controller) via Bluetooth or USB
   - Use apps like joystickMonitor to verify laptop can read controller
   - Terminal should display: "please calibrate joysticks"

4. **Power On Robot**
   - Connect 3S LiPo battery (voltage must be >11V)
   - Brushless motor controllers: red light → green light after ~3 seconds
   - Motors should twitch briefly

5. **Release Emergency Stop**
   - Twist to unpress the estop button if pressed
   - Front arm servos should move to default positions

6. **Zero Positions**
   - Press "circle" button on gamepad to zero and reset servo/encoder positions

7. **Verify Telemetry**
   - Check terminal output for telemetry data (printed around line 590 in brain2 code):
     - `v`: voltage
     - `d0-d4`: servo positions
     - `pA, pB`: BLDC powers
     - `mA, mB`: BLDC positions
     - `eA, eB`: leg encoder positions
     - `t`: time
     - `I`: current
   - If leg encoder positions are negative: manually retract both leg axes, then press circle button again

8. **Ready for Operation**
   - **Short forward jump:** Top of D-pad
   - **Long jump:** Left bumper + D-pad up
   - **Post-jump positioning:** Press 'B' and 'A' to move front arms for next jump setup

### Calibration Procedure

**When to Calibrate:**
- Leg encoders lose position (e.g., after disassembly)
- Motors continue running at endstop or stop too early during jumps
- Note: Issues may also be due to twisted string wound incorrectly or encoder failure

**Calibration Steps:**
1. Update encoder values in `poslib` dictionary in `station.py` (line 169)
2. Run `station.py`
3. Manually move leg links in unlatched mode
4. Record `eA` and `eB` values from terminal at three positions:
   - **"endstop":** Leg fully extended
   - **"retract":** Retracted but latch not reset
   - **"reset":** Fully retracted to reset latch
5. Divide terminal values by 100 and update `poslib`
6. Note: "trigger" entry is currently unused

## Research Publications

- **ICRA 2025 Paper:** "Tree-Jumping Robot" - [arXiv:2409.09203](https://arxiv.org/abs/2409.09203)
- **ICRA 2024 Presentation:** Force sensor jump analysis data
- **Data:** Experimental datasets available in `pintomech/ICRA24/` and `pintomech/ICRA25/`


## Open-source
Publically available repositories:
- https://github.com/qwertpas/pintocode
- https://github.com/qwertpas/pintomech
- Documentation: https://pintobotics.substack.com/p/documentation

This project is licensed under the MIT License - see the [LICENSE](pintocode/LICENSE) file for details.

---

**Contact:**
- Christopher Xu: christopher.y.x@gmail.com
- Jack Yan: jack_march_2019@outlook.com