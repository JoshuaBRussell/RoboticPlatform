from OpenGL.raw.GL.VERSION.GL_1_0 import GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, GL_PROJECTION, glOrtho
from PyQt5 import QtOpenGL
import OpenGL.GL as gl        # python wrapping of OpenGL
from OpenGL import GLU        # OpenGL Utility Library, extends OpenGL functionality

import math


class CoP_GL(QtOpenGL.QGLWidget):

    def __init__(self, parent):
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.cop_x = 0.1
        self.cop_y = 0.0

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

        self.update_cop_location()


    def update_cop_location(self):

        #Draw a a ring around the origin
        self.draw_circle(0.0, 0.0, 0.5, 255.0, 255.0, 255.0)
        self.draw_circle(0.0, 0.0, 0.45, 0.0, 0.0, 0.0)
        
        
        if math.sqrt(self.cop_x**2 + self.cop_y**2) < 0.45:
            self.draw_circle(self.cop_x, self.cop_y, 0.1, 0.0, 255.0, 0.0)
        else:
            self.draw_circle(self.cop_x, self.cop_x, 0.1, 255.0, 50.0, 0.0)
        


    def draw_circle(self, posx, posy, radius, r, g, b):
        gl.glColor(r/255.0, g/255.0, b/255.0) # set orange color
        sides = 200            
        gl.glLineWidth(5.0)
        gl.glBegin(gl.GL_POLYGON)    
        for i in range(200):    
            cosine= radius * math.cos(i*2*math.pi/sides) + posx 
            sine  = radius * math.sin(i*2*math.pi/sides) + posy 
            gl.glVertex2f(cosine,sine)
        gl.glEnd()

        
        