from typing import Any, Callable
from enum import Enum

from PyQt6.QtCore import Qt, QThread, pyqtSignal, QSize
from PyQt6.QtGui import (
    QFileSystemModel,
    QIcon,
    QMovie,
    QTextCursor
)
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QPushButton,
    QStyleFactory,
    QTextEdit,
    QVBoxLayout,
    QTreeView,
    QLabel,
    QHeaderView,
    QLineEdit,
    QProgressBar
)
import qdarkstyle


class Status(Enum):
    IDLE = "Idle"
    RUNNING = "Running"
    STOPPING = "Stopping"


class SetupThread(QThread):
    def __init__(self, parent, setup_func):
        super().__init__(parent)
        self.setup_func = setup_func

    def run(self):
        self.setup_func()


class MainDialog(QDialog):
    log_signal = pyqtSignal(str, dict)
    log_raw_signal = pyqtSignal(str)
    set_value_progress_bar1_signal = pyqtSignal(int, arguments=['value'])
    set_max_progress_bar1_signal = pyqtSignal(int, arguments=['value'])
    set_format_progress_bar1_signal = pyqtSignal(str, arguments=['format'])
    set_value_progress_bar2_signal = pyqtSignal(int, arguments=['value'])
    set_max_progress_bar2_signal = pyqtSignal(int, arguments=['value'])
    set_format_progress_bar2_signal = pyqtSignal(str, arguments=['format'])
    updateTreeView_signal = pyqtSignal()
    updateTreeView_done_signal = pyqtSignal()

    def __init__(self, app: QApplication):
        super(MainDialog, self).__init__(None)
        self.app = app

        self.setup_ui()

        self.set_status(Status.RUNNING)

        # Instantiate the SetupThread class and start the thread
        from processing import setup
        self.setup_thread = SetupThread(self, lambda: setup(self))
        self.setup_thread.finished.connect(self.on_setup_thread_finished)
        self.setup_thread.start()

        # Connect signals to slots
        self.log_signal.connect(self.log)
        self.log_raw_signal.connect(lambda msg: self.log(msg, {}))
        self.set_value_progress_bar1_signal.connect(self.set_value_progress_bar1)
        self.set_max_progress_bar1_signal.connect(self.set_max_progress_bar1)
        self.set_format_progress_bar1_signal.connect(self.set_format_progress_bar1)
        self.set_value_progress_bar2_signal.connect(self.set_value_progress_bar2)
        self.set_max_progress_bar2_signal.connect(self.set_max_progress_bar2)
        self.set_format_progress_bar2_signal.connect(self.set_format_progress_bar2)
        self.updateTreeView_signal.connect(self.updateTreeView)

    def setup_ui(self):
        self.setup_size()
        self.setup_icon()
        self.setup_window_flags()
        self.setup_palette()
        self.setup_status_bar()
        self.setup_main_layout()
        self.setWindowTitle("Genøme Project")
        self.changeStyle('Fusion')

    def setup_size(self):
        # Set size to 3/4 of the width of the screen and 3/4 of its height
        screen = QApplication.primaryScreen()
        size = screen.availableSize()
        self.resize(size.width() * 3 // 4, size.height() * 3 // 4)

    def setup_icon(self):
        self.setWindowIcon(QIcon('assets/icon.png'))
        self.setWindowIconText("Genøme Project")

    def setup_window_flags(self):
        # ensure that there is a minimize and maximize button
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.WindowMinMaxButtonsHint)

    def setup_palette(self):
        self.originalPalette = QApplication.palette()

    def setup_status_bar(self):
        # Create QLabel for the animated GIF
        self.runningGifLabel = QLabel()
        self.runningGif = QMovie("assets/running.gif")  # Load the GIF from "assets/running.gif"
        self.runningGifLabel.hide()

        self.statusLabel = QLabel("Status: <span style=\"color: gray; font-weight: bold;\">Idle</span>")

        # Get the height of the text and scale the GIF accordingly
        text_height = self.statusLabel.sizeHint().height()
        self.runningGif.setScaledSize(QSize(text_height, text_height))

        self.runningGifLabel.setMovie(self.runningGif)
        self.runningGif.start()

        label = QLabel("Created by Maxime DROUHIN, Orfeú MOURET, Maty Laye SOUARE and Maxime MAIRE.")

        self.statusBar = QHBoxLayout()
        self.statusBar.addWidget(self.statusLabel)
        self.statusBar.addWidget(self.runningGifLabel)
        self.statusBar.addStretch(1)
        self.statusBar.addWidget(label)

    def setup_main_layout(self):
        self.createTopLayout()
        self.createLeftGroupBox()
        self.createTopRightGroupBox()
        self.createProgressBar1GroupBox()
        self.createProgressBar2GroupBox()

        mainLayout = QGridLayout()
        mainLayout.addLayout(self.topLayout, 0, 0, 1, 2)
        leftLayout = QGridLayout()
        leftLayout.addWidget(self.LeftGroupBox, 0, 0)
        rightLayout = QGridLayout()
        rightLayout.addWidget(self.topRightGroupBox, 0, 0)
        rightLayout.addWidget(self.progressBar1GroupBox, 1, 0)
        rightLayout.addWidget(self.progressBar2GroupBox, 2, 0)
        mainLayout.addLayout(leftLayout, 1, 0)
        mainLayout.addLayout(rightLayout, 1, 1)
        mainLayout.addLayout(self.statusBar, 2, 0, 1, 2)
        mainLayout.setColumnStretch(0, 2)
        mainLayout.setColumnStretch(1, 3)

        self.setLayout(mainLayout)

        self.setWindowTitle("Genøme Project")
        self.changeStyle('Fusion')

    def on_setup_thread_finished(self):
        self.set_status(Status.IDLE)
        self.LeftGroupBox.setDisabled(False)

    def changeStyle(self, styleName: str):
        # https://chat.openai.com/c/8ab9e243-8e0c-4a9e-b9c8-fd53f0dcd61c
        QApplication.setStyle(QStyleFactory.create(styleName))
        QApplication.setPalette(QApplication.style().standardPalette())

        # Use QDarkStyle theme
        dark_style = qdarkstyle.load_stylesheet(qt_api='pyqt6')
        self.app.setStyleSheet(dark_style)

    def createTopLayout(self):
        self.topLayout = QHBoxLayout()

    def createTreeView(self):
        self.organisms_selection_label = QLabel("Organisms selection:")

        # Create an empty QTreeView
        self.treeView = QTreeView()

    def updateTreeView(self):
        # Create a QFileSystemModel and set its root directory to the current directory
        self.fileSystemModel = QFileSystemModel()
        self.fileSystemModel.setRootPath(".")

        # Set the model to the file system model
        self.treeView.setModel(self.fileSystemModel)

        # Hide the columns showing file sizes and file types and set the root index to the directory above the current directory
        self.treeView.hideColumn(1)
        self.treeView.hideColumn(2)
        self.treeView.setRootIndex(self.fileSystemModel.index("./Results"))
        # Make the first column take up all the available space
        self.treeView.header().setStretchLastSection(False)
        self.treeView.header().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)

        # Emit a signal when the updateTreeView is done
        self.updateTreeView_done_signal.emit()

    def createLeftGroupBox(self):
        self.LeftGroupBox = QGroupBox("Options")
        self.LeftGroupBox.setDisabled(True)

        self.createTreeView()

        functional_regions_label = QLabel("Functional regions selection:")
        self.checkBox1 = QCheckBox("CDS")
        self.checkBox2 = QCheckBox("centromere")
        # self.checkBox2.setDisabled(True)
        self.checkBox3 = QCheckBox("intron")
        self.checkBox4 = QCheckBox("mobile_element")
        # self.checkBox4.setDisabled(True)
        self.checkBox5 = QCheckBox("ncRNA")
        # self.checkBox5.setDisabled(True)
        self.checkBox6 = QCheckBox("rRNA")
        # self.checkBox6.setDisabled(True)
        self.checkBox7 = QCheckBox("telomere")
        # self.checkBox7.setDisabled(True)
        self.checkBox8 = QCheckBox("tRNA")
        # self.checkBox8.setDisabled(True)
        self.checkBox9 = QCheckBox("3'UTR")
        # self.checkBox9.setDisabled(True)
        self.checkBox10 = QCheckBox("5'UTR")
        # self.checkBox10.setDisabled(True)

        # Set tooltips for each checkbox
        self.checkBox1.setToolTip("coding sequence")
        self.checkBox2.setToolTip("chromosomal region")
        self.checkBox3.setToolTip("non-coding region within a gene")
        self.checkBox4.setToolTip("DNA sequence that can move within the genome")
        self.checkBox5.setToolTip("non-coding RNA molecule")
        self.checkBox6.setToolTip("ribosomal RNA molecule")
        self.checkBox7.setToolTip("chromosomal end")
        self.checkBox8.setToolTip("RNA molecule that carries amino acids to the ribosome")
        self.checkBox9.setToolTip("3' untranslated region: non-coding region downstream of the coding sequence")
        self.checkBox10.setToolTip("5' untranslated region: non-coding region upstream of the coding sequence")

        self.functional_regions = [
            self.checkBox1,
            self.checkBox2,
            self.checkBox3,
            self.checkBox4,
            self.checkBox5,
            self.checkBox6,
            self.checkBox7,
            self.checkBox8,
            self.checkBox9,
            self.checkBox10
        ]


        # Create a one-line text input area
        textInput_label = QLabel("Supplementary functional region:")
        self.textInput = QLineEdit()
        self.textInput.setPlaceholderText("ex: preRNA, D-loop, etc.")
        # self.textInput.setDisabled(True)


        # Create the main button
        self.mainButton = QPushButton("Start")
        self.mainButton.setDefault(True)
        self.mainButton.clicked.connect(self.on_main_button_pressed)


        layout = QGridLayout()

        layout.addWidget(self.organisms_selection_label, 0, 0, 1, 2)
        layout.addWidget(self.treeView, 1, 0, 1, 2)

        layout.addWidget(functional_regions_label, 2, 0, 1, 2)
        layout.addWidget(self.checkBox1, 3, 0)
        layout.addWidget(self.checkBox2, 4, 0)
        layout.addWidget(self.checkBox3, 5, 0)
        layout.addWidget(self.checkBox4, 6, 0)
        layout.addWidget(self.checkBox5, 7, 0)
        layout.addWidget(self.checkBox6, 3, 1)
        layout.addWidget(self.checkBox7, 4, 1)
        layout.addWidget(self.checkBox8, 5, 1)
        layout.addWidget(self.checkBox9, 6, 1)
        layout.addWidget(self.checkBox10, 7, 1)

        layout.addWidget(textInput_label, 8, 0, 1, 2)
        layout.addWidget(self.textInput, 9, 0, 1, 2)

        layout.addWidget(self.mainButton, 10, 0, 1, 2)

        # layout.addStretch(1)
        self.LeftGroupBox.setLayout(layout)

    def createTopRightGroupBox(self):
        self.topRightGroupBox = QGroupBox("Logs")

        self.textEdit = QTextEdit()
        self.textEdit.setReadOnly(True)

        layout = QVBoxLayout()
        layout.addWidget(self.textEdit)
        self.topRightGroupBox.setLayout(layout)

    def createProgressBar1GroupBox(self):
        self.progressBar1GroupBox = QGroupBox("Total progress (number of organisms treated)")

        self.progressBar1 = QProgressBar()
        self.progressBar1.setFormat("%v out of ??")

        layout = QVBoxLayout()
        layout.addWidget(self.progressBar1)
        self.progressBar1GroupBox.setLayout(layout)

    def createProgressBar2GroupBox(self):
        self.progressBar2GroupBox = QGroupBox("Current organism progress (percentage of NC files treated)")

        self.progressBar2 = QProgressBar()
        self.progressBar2.setFormat("%v%")

        layout = QVBoxLayout()
        layout.addWidget(self.progressBar2)
        self.progressBar2GroupBox.setLayout(layout)

    def set_value_progress_bar1(self, value):
        self.progressBar1.setValue(value)
    def set_max_progress_bar1(self, value):
        self.progressBar1.setMaximum(value)
    def set_format_progress_bar1(self, format):
        self.progressBar1.setFormat(format)

    def set_value_progress_bar2(self, value):
        self.progressBar2.setValue(value)
    def set_max_progress_bar2(self, value):
        self.progressBar2.setMaximum(value)
    def set_format_progress_bar2(self, format):
        self.progressBar2.setFormat(format)

    def set_status(self, status: Status):
        if status == Status.IDLE:
            self.statusLabel.setText("Status: <span style=\"color: gray; font-weight: bold;\">Idle</span>")
            self.runningGifLabel.hide()
        elif status == Status.RUNNING:
            self.statusLabel.setText("Status: <span style=\"color: green; font-weight: bold;\">Running…</span>")
            self.runningGifLabel.show()
        elif status == Status.STOPPING:
            self.statusLabel.setText("Status: <span style=\"color: orange; font-weight: bold;\">Stopping…</span>")
            self.runningGifLabel.show()

    def log(self, string: str, kwargs: dict):
        """
        Log a string to the textEdit widget.
        Here is a list of the available colors: https://htmlcolorcodes.com/color-names/
        """
        color = kwargs.get("color", "")
        bold = kwargs.get("bold", False)
        newline = kwargs.get("newline", True)

        # Move the text cursor to the end of the text
        cursor = self.textEdit.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.textEdit.setTextCursor(cursor)

        # Append the new string of text to the current text
        self.textEdit.textCursor().insertHtml(
            f"<span style=\"font-family: consolas; margin-bottom: 300%; color: {color}; font-weight: {'bold' if bold else ''};\">{string}</span>{'<br>' if newline else ''}"
        )

        # Scroll to the bottom of the textEdit widget
        scrollBar = self.textEdit.verticalScrollBar()
        scrollBar.setValue(scrollBar.maximum())

    def on_main_button_pressed(self):
        from processing import on_start_button_pressed, on_stop_button_pressed
        if self.mainButton.text() == "Start":
            self.mainButton.setText("Stop")
            self.set_status(Status.RUNNING)

            self.processing_thread = ProcessingThread(self, lambda: on_start_button_pressed(self), lambda: on_stop_button_pressed(self))
            self.processing_thread.finished.connect(self.on_processing_thread_finished)
            self.processing_thread.start()
        else:
            self.mainButton.setEnabled(False)
            self.set_status(Status.STOPPING)
            self.processing_thread.stop()

    def on_processing_thread_finished(self):
        self.mainButton.setText("Start")
        self.set_status(Status.IDLE)
        self.mainButton.setEnabled(True)


class ProcessingThread(QThread):
    def __init__(self, parent: MainDialog, on_start_func: Callable[[], None], on_stop_func: Callable[[], None]):
        super().__init__(parent)
        self.on_start_func = on_start_func
        self.on_stop_func = on_stop_func
        self.stop_requested = False

    def run(self):
        self.on_start_func()
        if self.stop_requested:
            self.on_stop_func()

    def stop(self):
        self.stop_requested = True
