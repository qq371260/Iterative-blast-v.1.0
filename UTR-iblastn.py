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
        self.root.title('UTR-Blast')
        self.root.geometry('500x400+400+200')

        # 变量接收
        self.db_var = ttkbs.StringVar()
        self.query_var = ttkbs.StringVar()
        self.orf_length_var = ttkbs.StringVar()
        self.parameters_var = ttkbs.StringVar()
        self.label_var = ttkbs.StringVar()
        self.num = 0

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

    # 布局设置
    def setting(self):
        self.Frame1 = tk.Frame(self.root, background='black', bg='blue')

        # 界面设置
        self.blank = ttkbs.Label(self.Frame1, width=4, border=10)
        self.Label1 = ttkbs.Label(self.Frame1, text='db file:')
        self.Label2 = ttkbs.Label(self.Frame1, text='query file:')
        self.Label3 = ttkbs.Label(self.Frame1,text='ORF length')
        self.Label4 = ttkbs.Label(self.Frame1, text='parameters:')
        self.Label_notice = ttkbs.Label(self.Frame1,font=('微软雅黑',12),textvariable=self.label_var)
        self.Entry_user = ttkbs.Entry(self.Frame1, textvariable=self.db_var)
        self.Entry_password = ttkbs.Entry(self.Frame1, textvariable=self.query_var)
        self.Entry_orf_length = ttkbs.Entry(self.Frame1,textvariable=self.orf_length_var)
        self.Entry_parameters = ttkbs.Entry(self.Frame1,textvariable=self.parameters_var)
        self.button_go = ttkbs.Button(self.Frame1, text='Go', bootstyle=PRIMARY, width=10,command=self.go)

        # 网格布局
        self.Frame1.grid()
        self.blank.grid(row=0, column=0, pady=25)
        self.Label1.grid(row=1, column=1, pady=4)
        self.Label2.grid(row=2, column=1, pady=4)
        self.Label3.grid(row=3, column=1, pady=4)
        self.Label4.grid(row=4,column=1,pady=4)
        self.Label_notice.grid(row=6, column=0,columnspan=4)
        self.Entry_user.grid(row=1, column=2, pady=3)
        self.Entry_password.grid(row=2, column=2, pady=4)
        self.Entry_orf_length.grid(row=3,column=2,pady=4)
        self.Entry_parameters.grid(row=4,column=2,pady=4)
        self.button_go.grid(row=5, column=2, pady=9)

    # Delete coding region nucleotide function
    def seek(self,filename):
        sequences = []
        descr = None
        # here is the path of multifalsta file
        with open(f"{filename}",encoding='utf-8') as file:
            line = file.readline()[:-1]  # always trim newline
            while line:
                if line[0] == '>':
                    if descr:  # any sequence found yet?
                        sequences.append((descr, seq))
                    descr = str(line[0:])
                    seq = ''  # start a new sequence
                else:
                    seq += line
                line = file.readline()[:-1]
            if descr:
                sequences.append((descr, seq))

        # 存储orf
        listOfOrf = list()
        # 存储剔除编码区后的核苷酸序列
        sequences_data = []
        for index, value in enumerate(sequences):  # looping over the fragments extracted
            global frames
            frames = [] # storing the six frame translation that it zould be extacted from the fragments
            dna = value[1]  # extract the fragment
            description = value[0] #extact the desciption even were not use it, just for learning purpose
            reverseCdna = [] # storing the reverse compliments
            # create the positive frames
            # split the frames into codons for better performance
            frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
            frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)])
            frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)])
            # reverse compliment of the fragment
            reverse = {"A": "T", "C": "G", "T": "A", "G": "C"}
            for i in range(len(dna)):
                reverseCdna.append(reverse[dna[-i - 1]]) if dna[-i - 1] in reverse.keys() else reverseCdna.append(dna[-i - 1])  # if any contamination found we keep it for further more check
            reverseCdna = ''.join(reverseCdna) # joining
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
                                listOfOrf.append(frames[i][start:stop+1])  # retrieve the orf
                                start = stop + 1  # avoiding multiple start codons
                                break
                    start += 1
            frame_orf = []
            # 存储待定位删除的数据
            orf_list_for_delete = []
            for i in listOfOrf:
                orf = ''.join(i)
                if len(i) > int(self.orf_length_var.get()):
                    orf_list_for_delete.append(orf)
                frame_cache = []
                for t in range(len(orf)):
                    frame_cache.append(reverse[orf[-t - 1]]) if orf[-t - 1] in reverse.keys() else frame_cache.append(orf[-t - 1])
                frame_cache = ''.join(frame_cache)
                if len(frame_cache) > int(self.orf_length_var.get())*3:
                    frame_orf.append(frame_cache)
            orf_list_for_delete.extend(frame_orf)
            # 遍历每一个核苷酸序列  变量dna
            # 还有名称 记录下来

            # 存储所有位置
            init_set = []
            # 遍历所有待删序列
            for j in orf_list_for_delete:
                index = re.compile(j)
                # 查找待删序列在核苷酸序列中的位置
                object = index.search(str(dna))
                if object:
                    init_set.append(object.span())
            total = []

            # 创建位置列表
            for index, i in enumerate(init_set):
                exec(f'self.list{index} = [x for x in range({i[0]},{i[1]})]')

            # 将所有位置载入列表total
            for t in range(len(init_set)):
                exec(f'total.extend(self.list{t})')
            section = []

            # 去重后加入到section中
            for i in range(len(total)):
                listx = [total[j] for j in range(i + 1, len(total))]
                if total[i] not in listx:
                    section.append(total[i])

            section.sort(reverse=True)
            # 根据section剔除相应位置
            for index in section:
                dna = dna[:index] + dna[index+1:]
            description = description
            sequences_data.append((description,dna))
        # 将剔除编码区后的核苷酸序列重新写入文件
        with open(f"{filename}", "w", encoding='utf-8') as file:
            for name, sequences in sequences_data:
                file.writelines(str(name) + "\n")
                file.writelines(str(sequences) + "\n")

    # cmd指令调整
    def blast(self,db_file,blast_file,c_db,parameters=''):
        self.ret = {}
        cmd1 = f'makeblastdb -in {db_file} -dbtype nucl -out {db_file}.blastdb'
        cmd2 = f'blastn -query {blast_file} -db {db_file}.blastdb -out {c_db} -outfmt 6 {parameters}'
        p = Popen(cmd1 + '&&' + cmd2, stdin=PIPE, stdout=PIPE,
                             stderr=PIPE, shell=True)
        p.wait()

    # 读取文件操作
    def read_file(self):
        # 读取B文件
        with open(self.c_db, 'r', encoding='utf-8') as fp:
            content = fp.read()
            self.Flag1 = content
            columns = re.findall(r'(^\w.*?)\t', content, re.MULTILINE)
            for column in columns:
                self.name.add(column)

        # 读取A文件
        self.Flag1 = os.path.getsize(f'./{self.blast_file}')
        with open(self.blast_file, 'r', encoding='utf-8') as fp:
            data = fp.read()
            name_pre = re.findall(r'^>(.*)', data, re.MULTILINE)
            result = re.findall(r'>.[^>]*', data, re.S)
            self.data_dict = {name_pre: result for name_pre, result in zip(name_pre, result)}

        if self.num >= 3:
            new_name, fa = os.path.splitext(str(self.blast_file))
            int_ = re.compile('\d')
            new_name = re.sub(int_, '', new_name)
            name = new_name + str(self.num-2) +  fa
            self.Flag2 = os.path.getsize(f'./{name}')

    # 修改文件操作
    def clean(self,file_c,file_a):
        # 创建新的B文件
        with open(file_c, 'w', encoding='utf-8') as fp:
            for i in self.name:
                if i in self.data_dict.keys():
                    fp.write(self.data_dict[i])
                    self.save_data[i] = self.data_dict[i]
                    del self.data_dict[i]

        # 创建新的A文件
        with open(file_a, 'w', encoding='utf-8') as fp:
            for i in self.data_dict.values():
                fp.write(i)

    # 文件名加数字函数
    def add_file_num(self,name):
        new_name, fa = os.path.splitext(name)
        int_ = re.compile('\d')
        new_name = re.sub(int_,"",new_name)
        name = new_name + str(self.num) +  fa
        return name

    # 制作数据库名称函数
    def db_name_maker(self):
        self.c_db = self.db_file
        c_db, fa = os.path.splitext(self.c_db)
        self.c_db = c_db + str(self.num) + '.blast'

    def c_blast_name(self):
        c_db, fa = os.path.splitext(self.c_db)
        self.c_db = c_db + str(self.num) +  fa

    # 替换R|Y|M|K|S|W|H|B|V|D特殊字符为N的函数
    def sub(self,file_name):
        with open(f"{file_name}", encoding='utf-8') as file:
            strings = file.read()

        file_name = self.add_file_num(file_name)
        with open(f"{file_name}", "w", encoding="utf-8") as file:
            strings = re.sub('R|Y|M|K|S|W|H|B|V|D', 'N', strings)
            file.write(strings)

    # 按键响应,流程控制函数
    def go(self):
        # 获取输入值
        self.db_file =self.db_var.get()
        self.blast_file = self.query_var.get()
        self.db_name_maker()
        self.parameters = self.parameters_var.get()

        # 开始流程
        try:
        # 替换所有R|Y|M|K|S|W|H|B|V|D
            self.sub(self.blast_file)
            self.sub(self.db_file)
            self.db_file = self.add_file_num(self.db_file)
            self.blast_file = self.add_file_num(self.blast_file)
            # 删除所有非编码区
            self.seek(self.blast_file)
            self.seek(self.db_file)
            if os.path.getsize(self.blast_file) and os.path.getsize(self.db_file):
                while  self.Flag1 != 0 and self.Flag2 != self.Flag1:
                    self.num += 1
                    # 判断有无参数，并执行cmd语句
                    if self.parameters:
                        self.blast( self.db_file, self.blast_file, self.c_db,parameters=self.parameters)
                    else:
                        self.blast( self.db_file, self.blast_file,self.c_db)
                    # 先读数据（便于后续对比）
                    self.read_file()
                    # 加数字
                    self.db_file = self.add_file_num(self.db_file)
                    self.blast_file = self.add_file_num(self.blast_file)
                    # 数据清洗，生成新的fa文件
                    self.clean(self.db_file,self.blast_file)

                    self.count += 1
                    # 执行完成
                # 创建B_total.txt文件
                with open('B_total.txt','w',encoding='utf-8') as fp:
                    for data in self.save_data.values():
                        regex = re.compile('^>(.*)')
                        fp.write(regex.search(data).group(1) + "\n")
                self.label_var.set('Blast Successfully')
                messagebox.showinfo('Success','Blast Successfully')
                # 创建A_iterations.txt文件
                with open('A_iterations.txt','w',encoding='utf-8') as fp:
                    fp.write(f'count:{self.count-2}')

                name, fa = os.path.splitext(self.blast_file)
                int_ = re.compile('\d')
                new_name = re.sub(int_, "", name)
                if self.num > 2:
                    name1 = new_name + str(self.num - 2) + fa
                else:
                    name1 = new_name +  fa
                new_name = new_name + fa

                # 将名字提取出来放入最后的文件中去
                name_list = []
                with open(name1,encoding='utf-8') as fp:
                    data = fp.readlines()
                    for i in data:
                        regex = re.compile('(>.*)')
                        name = regex.search(i)
                        if name:
                            name  = name.group(1)
                            name_list.append(name)
                with open(new_name, 'r', encoding='utf-8') as fp:
                    data_list = fp.readlines()
                    list1 = []
                    list2 = []
                    for line in data_list:
                        if line[0] == '>':
                            name = line[:-1]
                            seq = ''
                            if  seq:
                               list1.append(name)
                               list2.append(seq)
                        else:
                            seq += line
                    if name and seq:
                        list1.append(name)
                        list2.append(seq)
                    self.data_dict = {name: seq for name, seq in zip(list1, list2)}

                with open(self.blast_file,'w',encoding='utf-8') as fp:
                    for name in name_list:
                        if name in self.data_dict.keys():
                            fp.write(name + '\n')
                            fp.write(self.data_dict[name])

            else:
                messagebox.showinfo('Notice',
                f"No nucleotide sequence greater than {self.orf_length_var.get()} in the file {self.db_file} or {self.blast_file}")
        except Exception as e:
            messagebox.showinfo('Error',e)
            self.label_var.set('')

if __name__ == '__main__':
    Blast = Blast()