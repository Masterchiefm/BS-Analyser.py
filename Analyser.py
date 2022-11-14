import os.path

import markdown
from Bio import SeqIO
import regex, re
from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
import pandas as pd
from io import BytesIO
import base64

from Bio import pairwise2
from matplotlib import pyplot as plt


class SangerBaseCall:
    def __init__(self, sanger_file):
        self.sanger_data = SeqIO.read(sanger_file, "abi")
        self.sanger_seq = str(self.sanger_data.seq).upper()

        self.labels = {}
        signal_length = self.sanger_data.annotations["abif_raw"]["PLOC1"][-1]
        base_locations = self.sanger_data.annotations["abif_raw"]["PLOC1"]
        self.annotation = {}
        self.location_base_annotation = {}

        # 列出每个碱基位置上分别是什么碱基和他对应的测序峰值位置。
        for i in range(len(self.sanger_seq)):
            self.annotation[i] = [self.sanger_seq[i], base_locations[i]]
            self.location_base_annotation[base_locations[i]] = self.sanger_seq[i]

        # 列出每个信号列中的信号值
        channels = {"G": "DATA9", "A": "DATA10", "T": "DATA11", "C": "DATA12"}  # GATC
        self.trace = {}
        for c in channels:
            self.trace[c] = self.sanger_data.annotations["abif_raw"][channels[c]]

        # 根据信号峰值填入碱基位置
        for i in range(signal_length):
            if i in base_locations:
                self.labels[i] = self.location_base_annotation[i]
            else:
                self.labels[i] = ""

        # 整体的噪音
        self.noise = {"G": [], "A": [], "T": [], "C": []}

    def locTartet(self, target_seq, max_mismatch=10):
        # Find the match site and location
        target_is_reversed = False
        target_seq = str(target_seq).upper()
        alignment_1 = pairwise2.align.localms(target_seq, self.sanger_seq, 2, -.1, -4, -2, one_alignment_only=True)
        target_seq_reversed = str(Seq(target_seq).reverse_complement())
        alignment_2 = pairwise2.align.localms(target_seq_reversed, self.sanger_seq, 2, -.1, -4, -2,
                                              one_alignment_only=True)
        try:
            alignment_1[0]
        except:
            alignment_1 = [[0, 0, 0]]
        try:
            alignment_2[0]
        except:
            alignment_2 = [[0, 0, 0]]

        if alignment_2[0][2] > alignment_1[0][2]:
            alignment = alignment_2
        else:
            alignment = alignment_1
        # print(alignment)

        if target_is_reversed:
            location_start = alignment[0][4]
            location_end = alignment[0][3]
            # match_seq = str(alignment[location_end][location_start]).upper()
        else:
            location_start = alignment[0][3]
            location_end = alignment[0][4]
            match_seq = str(alignment[0][1][location_start:location_end]).upper()

        return (match_seq, location_start, location_end)


class BSReport:
    def __init__(self, sub_result_list):
        # 如果碱基发生变化，说明没有发生甲基化。0
        # 如果碱基没变，说明有甲基化，记为1
        # 将碱基转换成数字矩阵
        self.ref = sub_result_list[0]
        self.seqs = sub_result_list[1:]
        self.C_loc = []
        self.CG_loc = []
        for i in range(len(self.ref)):
            if self.ref[i] == "C":
                self.C_loc.append(i)
            if i > 0:
                if self.ref[i] == "G" and self.ref[i - 1] == "C":
                    self.CG_loc.append(i - 1)

    def getmCG(self):
        # print(self.seq)
        CG_peak = {}
        raw_data = []
        for loc in self.CG_loc:
            CG_peak[loc] = 0

        for seq in self.seqs:
            data = []
            for loc in self.CG_loc:
                if seq[loc] == self.ref[loc]:
                    CG_peak[loc] = CG_peak[loc] + 1
                    data.append(1)
                else:
                    data.append(0)
            raw_data.append(data)

        self.mCG_peak = CG_peak
        return (self.mCG_peak, raw_data)

    def getmC(self):
        # print(self.seq)
        C_peak = {}
        raw_data = []
        for loc in self.C_loc:
            C_peak[loc] = 0

        for seq in self.seqs:
            data = []
            for loc in self.C_loc:
                if seq[loc] == self.ref[loc]:
                    C_peak[loc] = C_peak[loc] + 1
                    data.append(1)
                else:
                    data.append(0)
            raw_data.append(data)

        self.mC_peak = C_peak
        return (self.mC_peak, raw_data)

    def plotPeakMap(self, peak_dic, title=""):
        plot_val = []
        for i in range(len(self.ref)):
            if i in peak_dic:
                plot_val.append(peak_dic[i])
            else:
                plot_val.append(0)


        plt.cla()
        plt.close("all")
        plt.rcParams['figure.figsize'] = (15, 4)
        colors = "black"
        fig, ax = plt.subplots()
        ax.plot(plot_val, color=colors[0], label='G')
        plt.xlabel("Reference sequence")
        plt.ylabel('Nums of methylation')
        plt.title(title)
        return plt

    def plotHeatMap(self, peak_dic, add_locs=[], exc_locs=[], title="",rows = []):
        plot_loc = list(peak_dic.keys())

        for i in add_locs:
            try:
                plot_loc.append(i-1)
            except:
                pass

        for i in exc_locs:
            try:
                plot_loc.remove(i-1)
            except:
                pass

        plot_loc.sort()

        data = []

        for seq in self.seqs:
            sub_data = []
            for loc in plot_loc:
                # print(seq[loc],self.ref[loc])
                if seq[loc] == self.ref[loc]:
                    if self.ref[loc] == "C":
                        sub_data.append(1)
                        # print(1)
                    else:
                        sub_data.append(0)
                else:
                    sub_data.append(0)
            sub_data.append(0)
            data.append(sub_data)
        # data_0 = lst1 = [0] * len(plot_loc)
        # data.append(data_0)
        data = np.array(data)

        loc_show = []
        for loc in plot_loc:
            loc_show.append(loc + 1)


        plt.cla()
        plt.close("all")
        plt.rcParams['figure.figsize'] = (15, 4)
        fig, ax = plt.subplots()
        fig.tight_layout()

        ax.imshow(data, cmap="binary")

        ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
        ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
        ax.grid(which="minor", color="gray", linestyle='-', linewidth=2)
        ax.set_xticks(np.arange(len(plot_loc)), labels=loc_show)
        ax.set_yticks(np.arange(len(rows)), labels=rows)
        plt.xlabel("CG positions")
        plt.title(title)

        return plt

    def plt2md(self, plt):
        """将plt转换成markdown图片"""
        buffer = BytesIO()
        plt.savefig(buffer, format='png')
        # 转换base64并以utf8格式输出
        pic_base64 = base64.b64encode(buffer.getvalue()).decode('utf8')
        md = """![](data:image/png;base64,""" + pic_base64 + """)"""
        return md

    def generateReport(self, ref_seq='', sanger_names=range(2), add_locs=[], exc_locs=[], report_name='test',
                       save_path=''):
        C_peak, C_data = self.getmC()
        CG_peak, CG_data = self.getmCG()

        C_peak_map = self.plotPeakMap(C_peak, title="Methalation state of C")
        C_peak_map_md = self.plt2md(C_peak_map)
        CG_peak_map = self.plotPeakMap(CG_peak, "Methalation state of CG")
        CG_peak_map_md = self.plt2md(CG_peak_map)

        CG_heat_map = self.plotHeatMap(CG_peak, add_locs=add_locs, exc_locs=exc_locs,
                                       title="Methalation status of each sequencing samples",rows=sanger_names)
        CG_heat_map_md = self.plt2md(CG_heat_map)

        C_peak_loc_show = []
        CG_peak_loc_show = []
        for i in list(C_peak.keys()):
            C_peak_loc_show.append(i+1)

        for i in list(CG_peak.keys()):
            CG_peak_loc_show.append(i+1)

        C_sheet = pd.DataFrame(C_data, index=sanger_names, columns=C_peak_loc_show)
        CG_sheet = pd.DataFrame(CG_data, index=sanger_names, columns=CG_peak_loc_show)

        sequences = []
        for seq in self.seqs:
            sequence = list(seq)
            sequences.append(sequence)
        sequence_sheet = pd.DataFrame(sequences, index=sanger_names, columns=range(1,len(ref_seq) + 1))
        sequence_sheet.loc["Reference"] = list(ref_seq)



        if save_path == "":
            pass

        elif "/" == save_path[-1]:
            pass
        else:
            save_path = save_path + "/"

        try:
            os.mkdir(save_path + "Reports/")
        except:
            pass
        try:
            os.mkdir(save_path + "Reports/Sheets/")
        except:
            pass

        C_sheet.to_excel(save_path + "Reports/Sheets/" + report_name + "_mC.xlsx")
        CG_sheet.to_excel(save_path + "Reports/Sheets/" + report_name + "_mCG.xlsx")
        sequence_sheet.to_excel(save_path + "Reports/Sheets/" + report_name + "_sequence.xlsx")

        # 网页报告内容
        title = "# BS-Seq Report of " + report_name + "\n"
        sep = "\n\n-------\n\n"
        description = (
                    "\nThis BS-seq report shows the methalation state of the sequence " + ref_seq[:10] + "..." + ref_seq[
                                                                                                                -10:] + "\n\n" +
                    "The program work on *ONLY ONE STRAND!*, which is the strand you input as reference. If the reference is CG, and the sequencing data is CG, too, this CG site will be mark as *Methalated*, or 1. \n\n" +

                    "这份BS-seq报告的参考序列是" + ref_seq[:10] + "..." + ref_seq[-10:] + "\n\n" +
                    "该程序仅根据参考序列所在的链进行分析，也就是说，请务必保证参考序列与测序结果的序列是几乎一致。在分析过程中，如果参考序列的CG在测序结果中也是CG，那么这是个甲基化位点，将会在结果表格中记为1\n")

        CG_heat_map_title = "## Methalation state of each sequencing sample"

        CG_sites = []
        for i in list(CG_peak.keys()):
            CG_sites.append(i+1)
        C_sites = []
        for i in list(C_peak.keys()):
            C_sites.append(i + 1)



        CG_heat_map_description = ("\nThis fig shows each sequencing samle's CG methalation state.\n\n" +
                                   "The ploted sites are " + str(CG_sites) + " and " + str(
                    add_locs) + ". \n\n These sites," + str(exc_locs) + ", did not plot in the figure as you request." +
                                   "\n\n该图展示了各个测序文件显示的甲基化CG位点," +
                                   "进行分析的位置有 " + str(CG_sites) + " 以及 " + str(
                    add_locs) + "。\n\n但是根据你设置的参数，" + str(exc_locs) + "这些点没有在图中画出。\n\n   \n\n")

        peak_map_title = "## Methalation peaks of C and CG"
        peak_map_description = (
                    "\nThe program calculated all the sequencing data and generated the methalation peak map\n\n" +
                    "通过分析汇总各个测序文件，本程序生成了参考序列上的CG和C甲基化峰图\n\n   \n\n")

        data_title = "## Data sheets"
        data_description = (
                    "\nThese excel files provide methalation infomations of all sequencing file. You can use them to generate pictures like this report, or use them for other purpose.\n\n" +
                    "这些文件提供了所有测序文件计算出的甲基化情况，你可以用这些数据生成与本报表一模一样的图，也可以拿去其它分析软件做更进一步的分析。这些数据保存在Reports/Sheets目录中。\n\n    \n\n")

        C_data_url = "\n\n[mC data sheet](Sheets/" + report_name + "_mC.xlsx" + ")"
        CG_data_url = "\n\n[mCG data sheet](Sheets/" + report_name + "_mCG.xlsx" + ")"
        sequence_url = "\n\n[Sequence info sheet](Sheets/" + report_name + "_sequence.xlsx" +")"





        md = [title, description, sep,
              CG_heat_map_title, CG_heat_map_description, CG_heat_map_md, sep,
              peak_map_title, peak_map_description, CG_peak_map_md, C_peak_map_md, sep,
              data_title, data_description, C_data_url, CG_data_url, sequence_url]

        h5 = markdown.markdown("".join(md), extensions=["tables"], encoding="utf-8")

        # print(h5)
        report_path = save_path + "Reports/" + report_name + ".html"
        with open(report_path, "w",encoding="utf-8") as f:
            f.write("""<meta http-equiv="Content-Type" content="text/html;charset=utf-8"/>\n\n""" + h5)
        return report_path












