from subprocess import Popen,PIPE
import os
import re
import tkinter as tk
from ttkbootstrap.constants import *
import ttkbootstrap as ttkbs
from tkinter import messagebox
class Blast(object):
    # 初始化函数
    def __init__(self):
        self.root = ttkbs.Window()
        self.root.title('Blastx')
        self.root.geometry('450x400+400+200')

        # 变量接收
        self.db_var = ttkbs.StringVar()
        self.query_var = ttkbs.StringVar()
        self.parameters_var = ttkbs.StringVar()
        self.label_var = ttkbs.StringVar()
        self.combobox_var = ttkbs.StringVar()
        self.combobox_var.set('')
        self.num = 1

        # 数据接收
        self.data_dict = {}
        self.save_data = {}
        self.name = set()
        self.new_C_file_list = []
        self.Flag = False
        self.Flag1 = True
        self.Flag2 = False
        self.count = 0


        self.setting()
        self.root.mainloop()

    # 查找ORF并翻译成氨基酸
    def orf_translate(self,file_path):
        self.sequences = []
        descr = None
        # here is the path of multifalsta file
        with open(r"{}".format(file_path), encoding='utf-8') as file:
            line = file.readline()[:-1]  # always trim newline
            seq = ''
            while line:
                if line[0] == '>':
                    if descr and len(seq) > 150:  # any sequence found yet?
                        self.sequences.append((descr, seq))
                    descr = str(line[1:].split('>'))
                    seq = ''  # start a new sequence
                else:
                    seq += line
                line = file.readline()[:-1]
            if len(seq) > 150:
                self.sequences.append((descr, seq))
        if not self.sequences:
            return
        self.listOfOrf = list()

        for index, value in enumerate(self.sequences):  # looping over the fragments extracted
            global frames
            frames = []  # storing the six frame translation that it zould be extacted from the fragments
            dna = value[1]  # extract the fragment
            description = value[0]  # extact the desciption even were not use it, just for learning purpose
            reverseCdna = []  # storing the reverse compliments
            # create the positive frames
            # split the frames into codons for better performance
            frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
            frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)])
            frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)])
            # reverse compliment of the fragment
            reverse = {"A": "T", "C": "G", "T": "A", "G": "C"}
            for i in range(len(dna)):
                reverseCdna.append(reverse[dna[-i - 1]]) if dna[-i - 1] in reverse.keys() else reverseCdna.append(
                    dna[-i - 1])  # if any contamination found we keep it for further more check
            reverseCdna = ''.join(reverseCdna)  # joining
            # create the negative frames
            frames.append([reverseCdna[i:i + 3] for i in range(0, len(reverseCdna), 3)])
            frames.append([reverseCdna[i:i + 3] for i in range(1, len(reverseCdna), 3)])
            frames.append([reverseCdna[i:i + 3] for i in range(2, len(reverseCdna), 3)])

            for i in range(0, len(frames), 1):  # looping all the frames
                start = 0
                while start < len(frames[i]):  # looping each frame for start and stop codons
                    if frames[i][start] == "ATG" or frames[i][start] == "TTG" or frames[i][start] == "CTG" or frames[i][start] == "GTG":
                        for stop in range(start + 1, len(frames[i]), 1):
                            if frames[i][stop] == "TAA" or frames[i][stop] == "TAG" or frames[i][stop] == "TGA":
                                self.listOfOrf.append(frames[i][start:stop])  # retrieve the orf
                                start = stop + 1  # avoiding multiple start codons
                                break
                    start += 1
        self.data_list = []
        table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
            'NNN': '_'
        }

        for j in self.listOfOrf:
            protein = ""
            for i in j:
                if "N" in i:
                    protein += table['NNN']
                else:
                    protein += table[i]
            self.data_list.append(protein)

        with open(r"{}".format(file_path), "w", encoding='utf-8') as file:
            for index, protein in enumerate(self.data_list):
                file.writelines(">" + str(index) + "\n")
                file.writelines(str(protein) + "\n")

    # 布局设置
    def setting(self):

        # 页面设置
        self.Frame1 = tk.Frame(self.root, background='black', bg='blue')

        # 界面设置
        self.blank = ttkbs.Label(self.Frame1, width=4, border=10)
        self.Label1 = ttkbs.Label(self.Frame1, text='db file:')
        self.Label2 = ttkbs.Label(self.Frame1, text='query file:')
        self.Label3 = ttkbs.Label(self.Frame1, text='parameters:')
        self.Label_notice = ttkbs.Label(self.Frame1,font=('微软雅黑',12),textvariable=self.label_var)


        self.Entry_user = ttkbs.Entry(self.Frame1, textvariable=self.db_var)
        self.Entry_password = ttkbs.Entry(self.Frame1, textvariable=self.query_var)
        self.Entry_parameters = ttkbs.Entry(self.Frame1,textvariable=self.parameters_var)

        self.button_go = ttkbs.Button(self.Frame1, text='Go', bootstyle=PRIMARY, width=10,command=self.go)

        # 网格布局
        self.Frame1.grid()
        self.blank.grid(row=0, column=0, pady=25)
        self.Label1.grid(row=1, column=1, pady=4)
        self.Label2.grid(row=2, column=1, pady=4)
        self.Label3.grid(row=3, column=1, pady=4)
        self.Label_notice.grid(row=5, column=0,columnspan=4)
        self.Entry_user.grid(row=1, column=2, pady=3)
        self.Entry_password.grid(row=2, column=2, pady=4)
        self.Entry_parameters.grid(row=3,column=2,pady=4)
        self.button_go.grid(row=4, column=2, pady=9)

    # cmd指令调整
    def blast(self,db_file,blast_file,c_db,parameters=''):
        self.ret = {}
        cmd1 = f'makeblastdb -in {db_file} -dbtype prot -out {db_file}.blastdb'
        cmd2 = f'blastx -query {blast_file} -db {db_file}.blastdb -out {c_db} -outfmt 6 {parameters}'
        p = Popen(cmd1 + '&&' + cmd2, stdin=PIPE, stdout=PIPE,
                             stderr=PIPE, shell=True)
        p.wait()

    # 读取文件操作
    def read_file(self):
        # 操作B库
        with open(self.c_db, 'r', encoding='utf-8') as fp:
            content = fp.read()
            self.Flag1 = content
            columns = re.findall(r'(^\w.*?)\t', content, re.MULTILINE)
            for column in columns:
                self.name.add(column)

        # 操作A库
        self.Flag1 = os.path.getsize(f'./{self.blast_file}')
        with open(self.blast_file, 'r', encoding='utf-8') as fp:
            data = fp.read()
            name_pre = re.findall(r'^>(.*)', data, re.MULTILINE)
            result = re.findall(r'>.[^>]*', data, re.S)
            self.data_dict = {name_pre: result for name_pre, result in zip(name_pre, result)}

        new_name, fa = os.path.splitext(str(self.blast_file))
        int_ = re.compile('\d')
        new_name = re.sub(int_,'', new_name)
        if self.num >= 3:
            name = new_name + str(self.num-2) +  fa
            self.Flag2 = os.path.getsize(f'./{name}')

    # 修改文件操作
    def clean(self,file_b,file_a):
        # 创建新的B库
        with open(file_b, 'w', encoding='utf-8') as fp:
            for i in self.name:
                if i in self.data_dict.keys():
                    fp.write(self.data_dict[i])
                    self.save_data[i] = self.data_dict[i]
                    del self.data_dict[i]
        # 查找ORF并翻译
        self.orf_translate(file_b)
        # 创建新的A库
        with open(file_a, 'w', encoding='utf-8') as fp:
            for i in self.data_dict.values():
                fp.write(i)

    # 文件名加数字函数
    def add_file_num(self,name):
        new_name, fa = os.path.splitext(name)
        int_ = re.compile('\d')
        new_name = re.sub(int_,"",new_name)
        name = new_name + str(self.num) + fa
        return name

    def db_name_maker(self):
        self.c_db = self.db_file
        c_db, fa = os.path.splitext(self.c_db)
        self.c_db = c_db + str(self.num) + '.blast'

    # 制作数据库名称函数
    def c_blast_name(self):
        c_db, fa = os.path.splitext(self.c_db)
        self.c_db = c_db + str(self.num) +  fa

    # 按键响应,流程控制函数
    def go(self):
        self.db_file =self.db_var.get()
        self.blast_file = self.query_var.get()
        self.db_name_maker()
        self.parameters = self.parameters_var.get()

        # 开始流程
        try:
            # 替换R|Y|M|K|S|W|H|B|V|D特殊字符为N的函数
            with open(f"{self.blast_file}", encoding='utf-8') as file:
                strings = file.read()
            with open(f"{self.blast_file}","w",encoding="utf-8") as file:
                strings = re.sub('R|Y|M|K|S|W|H|B|V|D','N',strings)
                file.write(strings)

            while  self.Flag1 != 0 and self.Flag2 != self.Flag1:
                # 判断有无参数，并执行cmd语句
                if self.parameters:
                    self.blast( self.db_file, self.blast_file, self.c_db,parameters=self.parameters)
                else:
                    self.blast( self.db_file, self.blast_file, self.c_db)
                # 先读数据（便于后续对比）
                self.read_file()
                # 加数字
                self.db_file = self.add_file_num(self.db_file)
                self.blast_file = self.add_file_num(self.blast_file)
                # 数据清洗，生成新的fa文件
                self.clean(self.db_file,self.blast_file)
                self.num += 1
                self.count += 1
                # 执行完成
            # 创建B_total.txt文件
            with open('B_total.txt','w',encoding='utf-8') as fp:
                for data in self.save_data.values():
                    fp.write(data)
            messagebox.showinfo('Success','Blast Successfully')
            # 创建A_iterations.txt文件
            with open('A_iterations.txt','w',encoding='utf-8') as fp:
                fp.write(f'count:{self.count-2}')
                a = 'A_iterations.txt'
        except Exception as e:
            messagebox.showinfo('Error',e)
            self.label_var.set('')

if __name__ == '__main__':
    Blast = Blast()
