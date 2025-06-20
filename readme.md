# Pinto: Tree-Jumping Squirrel Robot

![Pinto Robot Image](readme_media/treejumpclipi2.gif)


Pinto performs a tree-jump (jumping from the ground to a tree) like a squirrel. The robot uses a novel mechanism for the rear leg that can switch between a stiff series elastic and parallel elastic mode, for both precision movement and a powerful jump using elastic energy storage.


**Features:**
- **Dual-mode rear leg:** Series and parallel-elastic modes for both precise control and powerful jumps
- **Twisted string:** High mechanical advantage and lightweight transmission ideal for jumping
- **Spring:** Near constant torque carbon fiber spring for efficient windup 
- **Spined grippers:** Fishhooks for grasping tree bark during perching

**Paper Abstract:** Arboreal environments challenge current robots but are deftly traversed by many familiar animal locomotors such as squirrels. We present a small, 450 g robot "Pinto" developed for tree-jumping, a behavior seen in squirrels but rarely in legged robots: jumping from the ground onto a vertical tree trunk. We develop a powerful and lightweight latched series-elastic actuator using a twisted string and carbon fiber springs. We consider the effects of scaling down conventional quadrupeds and experimentally show how storing energy in a parallel-elastic fashion using a latch increases jump energy compared to series-elastic or springless strategies. By switching between series and parallel-elastic modes with our latched 5-bar leg mechanism, Pinto executes energetic jumps as well as maintains continuous control during shorter bounding motions. We also develop sprung 2-DoF arms equipped with spined grippers to grasp tree bark for high-speed perching following a jump.

## Repository Structure

### `pintocode/`
Contains firmware and control code for various microcontroller platforms and iterations:

#### **Version 1: Teensy 4.0 with Thumbstick Control**
- **`teensy/thumbstick*/`**: Early single-leg prototype controlled by thumbstick
- **Platform:** Teensy 4.0 microcontroller
- **Files:** `.ino` files for Arduino IDE

![Version 1 Image](readme_media/hopbox.gif)

#### **Version 2: ESP8266 Wireless Control System**
- **`teensy/joystick*/`**: Teensy 4.0 code for custom 8-axis joystick (four 2-axis thumbsticks)
- **`esp8266/espnow_send/`**: ESP8266 transceiver code for control station
- **`esp8266/espnow_recv/`**: ESP8266 receiver code for robot
- **Communication:** ESP-NOW wireless protocol
- **Control:** Human input via 8-axis joystick mapped to 4 robot limbs

<!-- ![Version 2 Image](readme_media/presentingv2.jpg) -->
<img src="readme_media/presentingv2.jpg" alt="Version 2 Image" width="480"/>


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


<p align="center">
<img src="readme_media/anatomy.png" alt="anatomy" width="800"/>
</p>

<p align="center">
<img src="readme_media/blockdiag.png" alt="v3 system diagram" width="800"/>
</p>

<p align="center">
<img src="readme_media/pintosetup.jpeg" alt="setup" width="600"/>
<img src="readme_media/estopinside.jpeg" alt="estop inside" width="240"/>
</p>




### `pintomech/`
Contains mechanical analysis, simulation, and experimental data:

#### **Mechanism Analysis:**
- **`SEA/`**: Series Elastic Actuator analysis and optimization python notebooks: see https://pintobotics.substack.com/p/extra-variable-mechanical-advantage
- **`spring_theory/`**: Theoretical spring mechanism analysis
- **`spring_experiments/`**: Physical spring testing and validation data
- **`string_experiments/`**: Twisted string actuator characterization: see https://pintobotics.substack.com/i/137821874/twisted-string-experiment
- **`lambda_linkage/`**: Lambda-linkage mechanism analysis for a vertically moving boom for testing a jumping leg: see https://pintobotics.substack.com/p/extra-linkage-optimization
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
- **`logs/`**: Experimental logs and data files

<p align="center">
<img src="readme_media/forcejump.png" alt="jump experiment" width="800"/>
</p>



## Hardware

- **Custom PCBs:** Available at [github.com/qwertpas/squirrelbrain](https://github.com/qwertpas/squirrelbrain)
- **CAD Files:** Available at [Onshape CAD](https://cad.onshape.com/documents/0c72dae6c9475dd41cab7700/w/2583ef1f629f03790729c107/e/f951a8c75300fed560aa6177)
- **Microcontrollers:** Teensy 4.0, ESP8266, ESP32S3 (Seeedstudio XIAO)

### **Motors:**
- **Brushless Motors:** 2x HGLRC SPECTER 2105.5 Brushless Motor 2650KV from Aliexpress
- **Servo Motors:** 5x Dynamixel XL330-M077-T

### **Sensors:**
- **Encoders:** 2x MA702 on custom PCB: https://github.com/qwertpas/o12encoder
- **IMU:** 1x LSM6DSV on custom PCB: https://github.com/qwertpas/o12imu
- **Force Sensors:** (not on the robot) 2x 5kg load cells

### **Power:**
- **Battery:** Tattu 3S 450mAh or 500mAh LiPo 
- **Power Management:** squirrelbrain PCB: https://github.com/qwertpas/squirrelbrain

### **Additional Components:**
- **Springs:** 0.5mm thick pultruded carbon fiber strips from Amazon
- **String:** UHMWPE rope (McMaster-Carr 4377N2)
- **Grippers:** Fishhooks from Amazon
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

## Publication

- **ICRA 2025:** "Pinto: A latched spring actuated robot for jumping and perching" by Christopher Y. Xu, Jack Yan, and Justin K. Yim [arXiv:2409.09203](https://arxiv.org/abs/2409.09203)



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
