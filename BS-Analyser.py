from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication, QTableWidgetItem, QFileDialog, QMessageBox
import pandas as pd
import os, requests
from Bio import SeqIO

from gui import Ui_MainWindow

from Analyser import SangerBaseCall, BSReport

def getLyric():
    try:
        url2 = 'https://v1.jinrishici.com/all'
        lyric = requests.get(url2, timeout=1).json()
        content = lyric['content']
        try:
            origin = lyric['origin']
        except:
            origin = "Unknown"
        try:
            author = lyric['author']
        except:
            author = "Unknown"


        output = content + "\n\t\t\t" + "——《" + origin + "》\t" + author
        return output
    except Exception as e:
        output = "BS-Seq Analyser" + "\n\t\t\t" + "——Written by M.Q. at ShanghaiTech University"
        return output

class MyMainWin(QMainWindow, Ui_MainWindow):
    def __init__(self, parent = None):
        super(MyMainWin, self).__init__(parent)

        self.setupUi(self)
        self.version = "1.1.2"

        # 表格标题读取
        col_count = self.tableWidget.columnCount()
        self.col_names = []
        self.col_name_locations = {}
        for i in range(col_count):
            name = self.tableWidget.horizontalHeaderItem(i).text()
            self.col_names.append(name)
            self.col_name_locations[name] = i

        print(self.col_names)

        self.tableWidget.selectColumn(0) #设选择第一页

        self.tableWidget.itemSelectionChanged.connect(self.showSelection)

        self.pushButton_del_lines.clicked.connect(self.delLine)
        self.pushButton_add_line.clicked.connect(self.addLine)
        self.pushButton_clear_table.clicked.connect(self.clearTable)
        self.pushButton_export_sheet.clicked.connect(self.exportSheet)
        self.pushButton_import_from_sheet.clicked.connect(self.importFromSheet)
        self.pushButton_open_folder.clicked.connect(self.openFolder)
        self.pushButton_open_report.clicked.connect(self.openReport)

        self.plainTextEdit_ref_rec_auto.dropped.connect(self.InputRef)
        self.plainTextEdit_sanger_rec_auto.dropped.connect(self.InputSanger)

        self.tableWidget.clicked.connect(self.disableAutoFill)

        self.pushButton_start.clicked.connect(self.annalyse)

    def annalyse(self):

        self.label_lyric.setText("❤❤❤❤❤❤   Analyzing... Please waite!   ❤❤❤❤❤❤")
        save_path = self.setSavePath()
        lyric = getLyric()
        if save_path:
            pass
        else:
            self.label_lyric.setText(lyric)
            return

        table = self.tableWidget
        # sheet = pd.DataFrame(columns=self.col_names)

        for row in range(table.rowCount()):
            try:
                name = table.item(row, 0).text().strip()
                ref_sequence = table.item(row, 1).text().strip()
                sanger_files = table.item(row, self.col_name_locations["Sequencing Files"]).text().split(",")
                sanger_names = []
            except:
                continue

            try:
                add_locs = []
                add_locs_text = table.item(row,self.col_name_locations["Additional locations"]).text()
                add_locs_text = add_locs_text.replace("，",",")
                for loc in add_locs_text.split(","):
                    add_locs.append(int(loc))
            except:
                add_locs = []

            try:
                exc_locs = []
                exc_locs_text = table.item(row,self.col_name_locations["Excluded locations"]).text()
                exc_locs_text = exc_locs_text.replace("，",",")
                for loc in exc_locs_text.split(","):
                    exc_locs.append(int(loc))
            except:
                exc_locs = []


            sub_result = [ref_sequence]
            for sanger_file in sanger_files:
                sanger_name = str(os.path.basename(sanger_file)).split(self.lineEdit_sep.text())[0]
                sanger_names.append(sanger_name)
                seq = SangerBaseCall(sanger_file)
                alignment = seq.locTartet(ref_sequence)[0]
                sub_result.append(alignment)

            report = BSReport(sub_result)
            report_path = report.generateReport(ref_seq=ref_sequence,sanger_names=sanger_names,add_locs=add_locs,exc_locs=exc_locs,report_name=name,save_path=save_path)
            self.tableWidget.setItem(row,self.col_name_locations["Report"],QTableWidgetItem(report_path))

        self.label_lyric.setText(lyric)
        QMessageBox.about(self,"Done","Analysis complete!")





    def disableAutoFill(self):
        self.checkBox_auto_fill_col.setChecked(False)


    def InputRef(self):
        lyric = getLyric()
        self.label_lyric.setText(lyric)

        file_list = self.plainTextEdit_ref_rec_auto.toPlainText().split("file:///")

        self.plainTextEdit_ref_rec_auto.clear()

        # msg = True
        for file in file_list:
            file_path = file.replace("\n", "")
            if os.path.isfile(file_path):
                if "dna" in file_path.lower()[-4:]:
                    file_name = os.path.basename(file_path)
                    sequence = str(SeqIO.read(file_path, "snapgene").seq).upper()
                    sample_name = file_name

                    exist = False
                    for row in range(self.tableWidget.rowCount()):
                        sample_exist = self.tableWidget.item(row, 0)
                        try:
                            if sample_name == sample_exist.text():
                                self.tableWidget.setItem(row, self.col_name_locations["Ref Sequence"],
                                                         QTableWidgetItem(sequence))
                                exist = True
                        except:
                            continue

                    if exist == False:
                        newRow = self.tableWidget.rowCount() + 1
                        self.tableWidget.setRowCount(newRow)
                        self.tableWidget.setItem(newRow - 1, self.col_name_locations["Ref Name"],
                                                 QTableWidgetItem(sample_name))
                        self.tableWidget.setItem(newRow - 1, self.col_name_locations["Ref Sequence"],
                                                 QTableWidgetItem(sequence))

        self.tableWidget.resizeColumnToContents(0)



    def InputSanger(self):
        file_list = self.plainTextEdit_sanger_rec_auto.toPlainText().split("file:///")
        self.plainTextEdit_sanger_rec_auto.clear()
        # msg = True

        try:
            sanger_file_list = self.tableWidget.item(self.tableWidget.currentRow(),self.col_name_locations["Sequencing Files"]).text()
        except:
            sanger_file_list = ''

        for file in file_list:
            file_path = file.replace("\n", "")
            if os.path.isfile(file_path):
                # file_name = os.path.basename(file_path)
                # sequence = str(SeqIO.read(file_path, "snapgene").seq).upper()
                if "ab1" in file_path.lower()[-4:]:
                    sanger_file_list = sanger_file_list + "," + file_path
        if sanger_file_list[0] == ",":
            sanger_file_list = sanger_file_list[1:]
        self.tableWidget.setItem(self.tableWidget.currentRow(), self.col_name_locations["Sequencing Files"], QTableWidgetItem(sanger_file_list))

    def openFolder(self):
        os.startfile(self.lineEdit_path.text())

    def openReport(self):
        try:
            selection = self.selected_rows
        except:
            return
        for i in selection:
            try:
                report_file = self.tableWidget.item(i, self.col_name_locations["Report"]).text()
                os.startfile(report_file)
            except Exception as e:
                QMessageBox.about(self,"ERROR","Report Not Found!\n" + str(e))


    def delLine(self):
        selection = self.selected_rows
        for i in selection:
            try:
                self.tableWidget.removeRow(selection[0])
            except Exception as e:
                print(e)


    def exportSheet(self,tem_save = False):
        lyric = getLyric()
        self.label_lyric.setText(lyric)

        if tem_save == True:
            file_path = self.lineEdit_path.text() + "/临时存储.xlsx"
        else:
            file_path, type = QFileDialog.getSaveFileName(self, "存", "", "excel(*.xlsx)")

        sheet = pd.DataFrame(columns = self.col_names)
        if file_path:
            pass
        else:
            return

        table = self.tableWidget

        for row in range(table.rowCount()):
            data = []
            for col in range(len(self.col_names)):
                try:
                    text = table.item(row, col).text()
                except:
                    text = ""
                data.append(text)
            sheet.loc[row] = data

        sheet.to_excel(file_path)

    def importFromSheet(self):
        file_path, type = QFileDialog.getOpenFileName(self, "导入", "", "excel(*.xlsx)")
        if file_path:
            pass
        else:
            return

        sheet = pd.read_excel(file_path, index_col=0)
        # count = self.tableWidget.rowCount()
        self.tableWidget.clearContents()
        self.tableWidget.setRowCount(len(sheet.index))

        for row in range(len(sheet.index)):
            data = sheet.iloc[row]
            for col in range(len(sheet.columns)):
                text = str(sheet.iloc[row][col])
                # if col == 0:
                #     if text == "nan":
                #         text = ""
                if text == "nan":
                    text = ""

                self.tableWidget.setItem(row, col, QTableWidgetItem(text))

    def setSavePath(self):
        save_path = QFileDialog.getExistingDirectory(self, "选路径")
        if save_path:
            pass
        else:
            return False
        self.lineEdit_path.setText(save_path)
        try:
            os.mkdir(save_path + "/" + "Reports")
        except:
            pass

        return save_path

    def addLine(self):
        current_count = self.tableWidget.rowCount()
        self.tableWidget.setRowCount(current_count + 1)
        new_row = current_count + 1
        for i in range(len(self.col_names)):
            self.tableWidget.setItem(new_row, i, QTableWidgetItem(""))

    def showSelection(self):
        selection = self.tableWidget.selectedIndexes()
        # self.currentSelectedIndex = selection
        # self.tabWidget.setCurrentIndex(0)
        rows = []
        names = []

        for i in selection:
            row = i.row()
            try:
                name = self.tableWidget.item(row, 0).text()
            except:
                name = ""

            if row in rows:
                pass
            else:
                rows.append(row)
                names.append(name)
        if len(rows) > 10:
            self.label_selection.setText("选中了" + str(len(rows)) + "个")
        else:
            self.label_selection.setText(str(names))
        self.selected_rows = rows

        col = self.tableWidget.currentColumn()
        if self.checkBox_auto_fill_col.isChecked():
            try:
                first_item = self.tableWidget.item(rows[0], col)
                first_data = first_item.text()
            except:
                first_data = ""

            for i in selection:
                row = i.row()
                self.tableWidget.setItem(row, col, QTableWidgetItem(first_data))


    def clearTable(self):
        self.tableWidget.clearContents()
        self.tableWidget.setRowCount(0)
        lyric = getLyric()
        self.label_lyric.setText(lyric)


if __name__ == "__main__":
    import sys

    # trans = QtCore.QTranslator()
    # trans.load("./en")


    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    # app.installTranslator(trans)
    win = MyMainWin()
    win.show()
    sys.exit(app.exec_())