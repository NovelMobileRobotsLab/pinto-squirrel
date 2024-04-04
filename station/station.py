import serial.tools.list_ports
import time
import pygame
import os
import pandas as pd
from timeit import default_timer as timer
from periodics import PeriodicSleeper
import numpy as np

os.environ['SDL_JOYSTICK_ALLOW_BACKGROUND_EVENTS'] = '1'
pygame.init()
screen = pygame.display.set_mode((50, 50))

port = "/dev/cu.usbmodem1101"
ser = None
# ports = serial.tools.list_ports.comports()

joysticks = {}
joy_data = {
    'leftx': 0,
    'lefty': 0,
    'rightx': 0,
    'righty': 0,
    'lefttrigger': 0,
    'righttrigger': 0,
    'A': 0,
    'B': 0,
    'X': 0,
    'Y': 0,
    '-': 0,
    'home': 0,
    '+': 0,
    'leftstickbutton': 0,
    'rightstickbutton': 0,
    'leftbumper': 0,
    'rightbumper': 0,
    'dpadup': 0,
    'dpaddown': 0,
    'dpadleft': 0,
    'dpadright': 0,
    'circle': 0,
}
axes_calibrated_dict = {
    'leftx+': False,
    'leftx-': False,
    'lefty+': False,
    'lefty-': False,
    'rightx+': False,
    'rightx-': False,
    'righty+': False,
    'righty-': False
}

global axes_calibrated
axes_calibrated = False

servos = [2048] * 5
motors = [0]*2
aux = 0

done = False
def main():
    global ser
    global servos, motors, aux

    traj_rightarm = pd.read_csv("data/backandforth.csv")

    send_handler = PeriodicSleeper(send_to_estop, 0.01)

    start_time = timer()
    while not done:
        elapsed = (timer() - start_time) * 1000

        while(ser is None):
            try:
                ser = serial.Serial(port, baudrate=115200, timeout=1, stopbits=serial.STOPBITS_TWO)
            except:
                ser = None
                print("plug in estop pls")
                time.sleep(0.1)
        
        handle_joysticks()

        speed = 3

        if(joy_data['X']):
            [servos[0], servos[1]] = interp(traj_rightarm, elapsed, ['dxl_pos[0]', 'dxl_pos[1]'], speed=speed) #maybe use period as parameter instead
            [servos[3], servos[2]] = np.array([4095, 4095]) - interp(traj_rightarm, elapsed, ['dxl_pos[0]', 'dxl_pos[1]'], speed=speed, phase_shift=0.5)
        elif(joy_data['Y']):
            [servos[0], servos[1]] = interp(traj_rightarm, elapsed, ['dxl_pos[0]', 'dxl_pos[1]'], speed=speed) #maybe use period as parameter instead
            [servos[3], servos[2]] = np.array([4095, 4095]) - interp(traj_rightarm, elapsed, ['dxl_pos[0]', 'dxl_pos[1]'], speed=speed, phase_shift=0.0)
        elif(joy_data['A']):
            servos = [2417, 2465, 1730, 1717, 2048] #lean forward
        elif(joy_data['B']):
            servos = [1848, 769, 3330, 2584, 2048] #lean forward
        elif(joy_data['dpaddown']):
            servos = [3114, 3172, 3649, 3199, 1107]
        elif(joy_data['dpadup']):
            servos = [3114, 3172, 3649, 3199, 1588]
        
        if(joy_data['rightbumper']):
            servos[4] += 50
            # servos[4] = 0
        elif(joy_data['leftbumper']):
            servos[4] -= 50
            # servos[4] = 4095
        # servos[0] = 2048 + joy_data['leftx']*2048 
        # servos[1] = 2048 + joy_data['lefty']*2048 
        # servos[2] = 2048 + joy_data['rightx']*2048 
        # servos[3] = 2048 + joy_data['righty']*2048 
            

        '''
        feedforward = ff_lut(enc - string_lut(motor_pos))
        motor = P*des_pos + feedforward
        '''

        motors[0] = joy_data['lefty']*511
        motors[1] = joy_data['righty']*511

        for i in range(5):
            servos[i] = min(max(int(servos[i]), 0), 4095)
        for i in range(2):
            motors[i] = min(max(int(motors[i]), -511), 511)
            if(abs(motors[i]) < 20):
                motors[i] = 0


        aux = 0
        if(joy_data['circle'] and joy_data['home']):
            aux |= (1 << 0) #set bit 0
        if(joy_data['circle']):
            aux |= (1 << 1) #set but 1


        recv_from_estop()
        
        time.sleep(0.01)

    send_handler.stop()



# print message cleanly to terminal, along with other station data
def print_message(message):
    global servos, motors, aux
    if(len(joysticks) == 0):
        print("gamepad disconnected")
    else:
        if(not axes_calibrated):
            print("please calibrate joysticks")
            print("move axes in a circle")
        else:
            joy_df = pd.DataFrame([joy_data]).round(2)
            print(joy_df[['leftx', 'lefty', 'rightx', 'righty']].to_string(index=False))

            print('s:', np.int_(servos))
            print('m:', np.int_(motors))
            # print('a:', format(aux, '05b'))
            print(f" a:{aux:05b}")
    print(f"aux: {aux}")
    print(message)
    print("_—————____—————____—————____—————\t")




def send_to_estop():
    global ser
    if(ser is None):
        return

    cmd_str = ""
    for i in range(5):
        cmd_str += f"s{i}:{servos[i]:05}\n"
    for i in range(2):
        cmd_str += f"m{i}:{motors[i]:05}\n"
    cmd_str += f"a:{aux:05b}"
    cmd_str += "#\t"

    # print(cmd_str)

    try:
        ser.write(cmd_str.encode())
    except Exception as e:
        print(e)
        ser = None
        print("Estop disconnected")


message = ""
delimiter = '\t'
def recv_from_estop():
    # print(ser.read_all().decode("utf-8", errors='ignore'), end=None)
    global ser
    global messagebuffer
    global message
    global axes_calibrated

    messagecount = 0

    try:
        if ser.in_waiting > 0:
            uarttext = ser.read_all().decode('utf-8', errors='ignore')
            
            ending=0
            while(uarttext):
                ending = uarttext.find(delimiter)
                if(ending == -1):
                    break

                message += uarttext[0:ending]

                print_message(message)
                messagebuffer = message

                messagecount += 1
                message = "" #clear message
                uarttext = uarttext[ending+len(delimiter):] #front of buffer used up

            message = uarttext #whatver is left over
    except Exception as e:
        print(e)
        return


def interp(traj_df, elapsed, labels, speed=1, phase_shift=0, dir=1):
    period = traj_df['elapsed'][len(traj_df)-1]/speed

    phase = elapsed/period + phase_shift
    progress = (phase - int(phase)) * len(traj_df) #goes from 0 to length of traj
    index = int(progress)
    index_frac = progress - index

    if(index + 1 < len(traj_df)):
        index_next = index + 1
    else:
        index_next = 0

    interpolated = []
    for label in labels:
        interpolated.append((1-index_frac)*traj_df[label][index] + (index_frac)*traj_df[label][index_next])
    return np.array(interpolated)
        


def handle_joysticks():
    global axes_calibrated
    global axes_calibrated_dict

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            global done
            done = True

        if event.type == pygame.JOYDEVICEADDED:
            joy = pygame.joystick.Joystick(event.device_index)
            joysticks[joy.get_instance_id()] = joy
            print(f"Joystick {joy.get_instance_id()} connected")
            axes_calibrated = False

        if event.type == pygame.JOYDEVICEREMOVED:
            del joysticks[event.instance_id]
            print(f"Joystick {event.instance_id} disconnected")

        for joystick in joysticks.values():
            joy_data['leftx'] = joystick.get_axis(0) 
            joy_data['lefty'] = -joystick.get_axis(1)
            joy_data['rightx'] = joystick.get_axis(2)
            joy_data['righty'] = -joystick.get_axis(3) # -1 to 1
            joy_data['lefttrigger'] = joystick.get_axis(4) # -1 to 1
            joy_data['righttrigger'] = joystick.get_axis(5) # -1 to 1
            joy_data['A'] = joystick.get_button(0)
            joy_data['B'] = joystick.get_button(1)
            joy_data['X'] = joystick.get_button(2)
            joy_data['Y'] = joystick.get_button(3)
            joy_data['-'] = joystick.get_button(4)
            joy_data['home'] = joystick.get_button(5) #make sure to disable Launchpad https://apple.stackexchange.com/questions/458669/how-do-i-disable-the-home-button-on-a-game-controller-in-macos
            joy_data['+'] = joystick.get_button(6)
            joy_data['leftstickbutton'] = joystick.get_button(7)
            joy_data['rightstickbutton'] = joystick.get_button(8)
            joy_data['leftbumper'] = joystick.get_button(9)
            joy_data['rightbumper'] = joystick.get_button(10)
            joy_data['dpadup'] = joystick.get_button(11)
            joy_data['dpaddown'] = joystick.get_button(12)
            joy_data['dpadleft'] = joystick.get_button(13)
            joy_data['dpadright'] = joystick.get_button(14)
            joy_data['circle'] = joystick.get_button(15)

            #seems like pygame needs to know the max axis value so need to move stick to max/min
            for axis_name in ['leftx', 'lefty', 'rightx', 'righty']:
                if(joy_data[axis_name] > 0.99):
                    axes_calibrated_dict[f"{axis_name}+"] = True
                elif(joy_data[axis_name] < -0.99):
                    axes_calibrated_dict[f"{axis_name}-"] = True
            
            #check if each axis calibrated, complete only all done and at home
            if(not axes_calibrated):
                axes_calibrated_try = True
                for calibration_label in axes_calibrated_dict:
                    if(axes_calibrated_dict[calibration_label] == False):
                        axes_calibrated_try = False
                        break
                for axis_name in ['leftx', 'lefty', 'rightx', 'righty']:
                    if(abs(joy_data[axis_name]) > 0.01):
                        axes_calibrated_try = False
                        break
                axes_calibrated = axes_calibrated_try
                if(axes_calibrated):
                    joystick.rumble(0.7, 0.0, 100) #indicate calibration done
                else:
                    for axis_name in ['leftx', 'lefty', 'rightx', 'righty']:
                        joy_data[axis_name] = 0

            break #assume only one joystick


if(__name__ == "__main__"):
    main()