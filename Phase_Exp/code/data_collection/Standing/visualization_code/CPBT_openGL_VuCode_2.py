'''
Name: max_lean_main.py
Description: Combine the GUI and the visualization together
Author: Vu Phan
Date: June 2th 2020
'''

import new_code # main visualization GUI is here
# pyOpenGL (graphics library) and PyQt (GUI library) packages 
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL import *
from PyQt5.QtOpenGL import *
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
# Misc. packages 
import sys,time
import numpy as np
import serial
import serial.tools.list_ports
import threading
import pandas as pd

# --- Constant (i.e., macros in the Arduino code) --- #
# Useful config. for the transmission frame
STARTFRAME = 0xaa # start byte
STOPFRAME = 0x55 # stop byte
RXFRAMESIZE = 18 # size of the receive frame = 18 (including start/stop)
RXFIRSTBYTE = 0 # index of the start byte
RXLASTBYTE = RXFRAMESIZE - 1 # index of the stop byte
# Transmission frame of 16 bytes (not including start/stop byte)
# Right platform
A0_BYTE1 = 1 # data from A0 pin
A0_BYTE2 = 2
A1_BYTE1 = 3 # data from A1 pin
A1_BYTE2 = 4
A2_BYTE1 = 5 # data from A2 pin
A2_BYTE2 = 6
A3_BYTE1 = 7 # data from A3 pin
A3_BYTE2 = 8
# Left platform
A4_BYTE1 = 9 # data from A4 pin
A4_BYTE2 = 10
A5_BYTE1 = 11 # data from A5 pin
A5_BYTE2 = 12
A6_BYTE1 = 13 # data from A6 pin
A6_BYTE2 = 14
A7_BYTE1 = 15 # data from A7 pin
A7_BYTE2 = 16

class MainWindow(QtWidgets.QMainWindow, new_code.Ui_MainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)
        self.showMaximized() # full screen setup

        # Init
        self.var_init() # initialize variables
        self.home() # initialize signaling functions

    ''' Variable initialization '''
    def var_init(self):
        # Serial config
        self.port = ''
        self.baudrate = 0
        self.ser = serial.Serial() # initialize serial communication
        
        # Experiment information
        self.subject_name = ''
        self.compliant_level = ''
        self.visual_condition = ''

        # Visualization parameters
        self.base = 0.54

        # --- Data frame settings --- #
        # Receive buffer
        self.RX_data_frame = []
        for i in range(RXFRAMESIZE):
            self.RX_data_frame.append(0) # init values
    
        # Pointer to read data from the frame
        self.frame_ptr = 0 # pointer of frame
        self.busy_flag = False # to enable/disable receiving data
        # false -> enable
        # true -> disable

    ''' Settings & Signal handling '''
    def home(self):
        ''' Visualization '''        

        # Set the viewport to window dimensions
        self.gl_widget.resizeGL(self.gl_widget.width(), self.gl_widget.height())
        # Timer to update visualization
        timer = QtCore.QTimer(self)
        timer.setInterval(20) # in ms
        # Updates Left Plot 
        timer.timeout.connect(self.gl_widget.updateGL)
        # Updates Right Plot
        timer.timeout.connect(self.gl_widget_2.updateGL)
        timer.start()

        # Timer to collect data
        self.data_timer = QtCore.QTimer(self)
        self.data_timer.setInterval(1) # 10ms
        self.data_timer.timeout.connect(self.update)


        if not self.ser.isOpen():
            # if open create variable of serial data 
            self.ser = serial.Serial(port='COM17',baudrate=115200)

            # Start the timer to collect data and update the graph
            self.data_timer.start()
        else:
            # Stop the timer before close the serial to avoid errors
            self.data_timer.stop()
            self.ser.close() # close the serial communication

        # Check if we have connected to the serial port
        if self.ser.isOpen():
            print("Connected")
        else:
            print("Disconnected")

    # Update function
    def update(self):

        # Check A2, A3, A6, A7
        
        if self.ser.isOpen():
            try:
                #--- Extract data from the transmision frame ---#
                if self.busy_flag == False:
                    temp = self.ser.read(1, timeout=1) # store in a temporary variable
                    self.RX_data_frame[self.frame_ptr] = ord(temp) # convert to hex value

                    # Check the start byte
                    if self.RX_data_frame[RXFIRSTBYTE] == STARTFRAME:
                        self.frame_ptr = self.frame_ptr + 1 # next position
                    else:
                        # do nothing if failed
                        pass

                    # Copy data after having a whole (CORRECT) frame
                    if (self.RX_data_frame[RXFIRSTBYTE] == STARTFRAME) and \
                       (self.RX_data_frame[RXLASTBYTE] == STOPFRAME):
                        self.busy_flag = True # disable receiving data

                        # --- Extract data frame --- #
                        # Right platform
                        y_r = (self.RX_data_frame[A2_BYTE1] << 8) + \
                              self.RX_data_frame[A2_BYTE2]
                        x_r = (self.RX_data_frame[A3_BYTE1] << 8) + \
                              self.RX_data_frame[A3_BYTE2]
                        # Left platform
                        y_l = (self.RX_data_frame[A6_BYTE1] << 8) + \
                              self.RX_data_frame[A6_BYTE2]
                        x_l = (self.RX_data_frame[A7_BYTE1] << 8) + \
                              self.RX_data_frame[A7_BYTE2]

                        # Clear the receive buffer
                        for i in range(RXFRAMESIZE):
                            self.RX_data_frame[i] = 0
                        self.busy_flag = False # enable receiving data again
                        self.frame_ptr = RXFIRSTBYTE # 0

                        # ---  Update the circle position --- #
                        # Move the left subject circle, scaled based on monitor 
                        self.gl_widget.subjectx = -(x_l/1000.0 - self.base)/0.8
                        self.gl_widget.subjecty = -(y_l/1000.0 - self.base)/0.8

                        # Move the right subject circle, scaled based on monitor 
                        self.gl_widget_2.subjectx = -(x_r/1000.0 - self.base)/0.8
                        self.gl_widget_2.subjecty = -(y_r/1000.0 - self.base)/0.8
                        
                    else:
                        pass # do nothing
                else:
                    pass # do nothing              
                                    
            except:
                #print(temp)
                print("Fault detected during communication !!!")
                self.frame_ptr = RXFIRSTBYTE # reset pointer if fault
                
        else: # self.ser.isOpen():
            # do nothing
            pass

            
# Standard run code for pyqt application
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    # Display GUI
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec_())
