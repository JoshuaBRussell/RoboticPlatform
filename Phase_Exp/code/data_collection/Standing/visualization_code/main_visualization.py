from PyQt5 import QtGui, QtCore, QtWidgets
import cop_and_weight_layout

class MainWindow(QtWidgets.QMainWindow, cop_and_weight_layout.Ui_MainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)
        self.showMaximized() # full screen setup

        # Set the viewport to window dimensions
        self.cop_widget.resizeGL(self.cop_widget.width(), self.cop_widget.height())
        self.weight_widget.resizeGL(self.weight_widget.width(), self.weight_widget.height())

        # Timer to update visualization
        timer = QtCore.QTimer(self)
        timer.setInterval(20) # in ms
        # Updates Plot 
        timer.timeout.connect(self.cop_widget.updateGL)
        timer.timeout.connect(self.weight_widget.updateGL)
        timer.start()



if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)

    #Display GUI
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec_())