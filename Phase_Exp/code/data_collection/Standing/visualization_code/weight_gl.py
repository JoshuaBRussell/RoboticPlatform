from OpenGL.raw.GL.VERSION.GL_1_0 import GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, GL_PROJECTION, glOrtho
from PyQt5 import QtOpenGL
import OpenGL.GL as gl        # python wrapping of OpenGL
from OpenGL import GLU        # OpenGL Utility Library, extends OpenGL functionality

import math


class Weight_GL(QtOpenGL.QGLWidget):

    def __init__(self, parent):
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.index = 0.0

        self.normalized_weight = 0.0
        self.x_pos = 0.0
        self.y_pos = 0.0

    def initializeGL(self):
        gl.glClearColor(0.0, 0.0, 0.0, 0.0)
        gl.glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    def resizeGL(self, w: int, h: int):
        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        gl.glOrtho(-1, 1, -1, 1, -50.0, 50.0)
        gl.glViewport(0, 0, w, h)
        print("Width = " + str(w))
        print("Height = " + str(h) + '\n---')


    def paintGL(self) -> None:
        gl.glClearColor(0.0, 0.0, 0.0, 0.0) # black screen
        gl.glClear(GL_COLOR_BUFFER_BIT) # clear color buffer and set color

        # self.index = (self.index + 0.02)%1
        # self.normalized_weight = self.index
        self.update_weight_bar()


    def update_weight_bar(self):

        # #Draw a a ring around the origin
        self.draw_bar(0.05, 1.0, 0.0, 0.45, 255.0, 255.0, 255.0)
        self.draw_bar(0.05, 1.0, 0.0, -0.45,255.0, 255.0, 255.0)


        if  abs(self.normalized_weight < 0.45):
            self.draw_bar(self.normalized_weight, 0.5, self.x_pos, self.y_pos,   0.0, 255.0, 0.0)
        else:
            self.draw_bar(self.normalized_weight, 0.5, self.x_pos, self.y_pos, 255.0,  50.0, 0.0)


    def draw_bar(self, height, width, x_pos, y_pos, r, g, b):
        gl.glColor(r/255.0, g/255.0, b/255.0)

        gl.glLineWidth(5.0)

        gl.glBegin(gl.GL_QUADS)
        gl.glVertex2f(-width/2 + x_pos,    0.0 + y_pos)
        gl.glVertex2f(-width/2 + x_pos, height + y_pos)
        gl.glVertex2f( width/2 + x_pos, height + y_pos)
        gl.glVertex2f( width/2 + x_pos,    0.0 + y_pos)
        gl.glEnd()
