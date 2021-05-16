'''
Name: max_lean_game.py
Description: Visualization for max leaning task implemented by OpenGL
Author: Vu Phan
Date: June 2th 2020
'''

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL import *
from PyQt5.QtOpenGL import *
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
from math import *
import sys,time
import numpy as np
from CPBT_openGL_VuCode_2 import *

''' Constant '''
# Initial position of the target    
INIT_TARGET_FRONT = 0.55 # front direction
INIT_TARGET_BACK = 0.35 # backward direction
TARGET_STEP = 0.02 # step to increase/decrease
TARGET_RADIUS = 0.1 # radius of the target circle
MARKED_TARGET_RADIUS = 0.02 # radius of the marked target circle
SUBJECT_POS_RADIUS = 0.08
AREA_OUTER_RADIUS = 0.5 # radius of the outer circular area surrounding the user's cursor
AREA_INNER_RADIUS = 0.48 # radius of the inner circular area surrounding the user's cursor 
# Scale
MAX_VALUE = 20 # Visualization show from -20 to 20 cm in all direction

''' Visualization '''
class glWidget(QGLWidget):
    global INIT_TARGET_FRONT, INIT_TARGET_BACK, TARGET_STEP, \
           TARGET_RADIUS, MARKED_TARGET_RADIUS, SUBJECT_POS_RADIUS, \
           MAX_VALUE, AREA_RADIUS

    def __init__(self, parent):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(400, 400)

        # Variable initialization
        self.var_init()

    def var_init(self):

        ''' Subject position '''
        self.subjectx = 0 # x-dir position of the subject <- controlled by GUI
        self.subjecty = 0 # y-dir position of the subject <- controlled by GUI

    def initializeGL(self):
        glClearColor(0.0, 0.0, 0.0, 0.0) # black screen
        glClear(GL_COLOR_BUFFER_BIT) # clear color buffer and set color

    def resizeGL(self, w, h):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity() # reset coordinate syste
        glOrtho(-w/h, w/h, -1, 1, -50.0, 50.0) # scale the matrix
            # (left, right, bottom, top, zNear, zFar)
        glViewport(0, 0, w, h) # set viewport to window dimensions
        print("Width = " + str(w))
        print("Height = " + str(h) + '\n---')

    def paintGL(self):
        glClearColor(0.0, 0.0, 0.0, 0.0) # black screen
        glClear(GL_COLOR_BUFFER_BIT) # clear color buffer and set color

        ''' Area '''
        self.gl_boundary()
        ''' Subject circle '''
        self.gl_subject_pos()

    ''' Subject position '''
    def gl_subject_pos(self):
        glColor(1.0, 165.0/255.0, 0.0) # set orange color
        self.draw_circle(self.subjectx, self.subjecty, SUBJECT_POS_RADIUS)

    ''' Boundary '''
    def gl_boundary(self):
        glColor(1,1,1)
        self.draw_outer_area(AREA_OUTER_RADIUS)
        glColor(0, 0, 0) # inner circle should match the background
        self.draw_inner_area(AREA_INNER_RADIUS)

    ''' Useful functions '''
    def draw_circle(self, posx, posy, radius):
        sides = 200            
        glLineWidth(5.0)
        glBegin(GL_POLYGON)    
        for i in range(200):    
            cosine= radius * cos(i*2*pi/sides) + posx 
            sine  = radius * sin(i*2*pi/sides) + posy 
            glVertex2f(cosine,sine)
        glEnd()

        ''' Useful functions '''
    def draw_outer_area(self,radius_outer):
        sides = 200            
        glBegin(GL_POLYGON)    
        for i in range(200):    
            cosine= radius_outer * cos(i*2*pi/sides) + 0
            sine  = radius_outer * sin(i*2*pi/sides) + 0 
            glVertex2f(cosine,sine)
        glEnd()

    def draw_inner_area(self, radius_inner):
        sides = 200            
        glBegin(GL_POLYGON)    
        for i in range(200):    
            cosine= radius_inner * cos(i*2*pi/sides) + 0
            sine  = radius_inner * sin(i*2*pi/sides) + 0 
            glVertex2f(cosine,sine)
        glEnd()