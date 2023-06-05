import sys
import os

from gui import MainDialog, QApplication


if __name__ == '__main__':
    # Set current directory to the directory of the script
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Starting the app. From now on, functions will be executed upon interaction
    # with the UI by the user.
    app = QApplication(sys.argv)
    main_dialog = MainDialog(app)
    main_dialog.show()
    sys.exit(app.exec())
