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
        self.root.title('BLASTN & TBLASTX')
        self.root.geometry('550x400+400+200')

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
        # 监听事件
        self.root.mainloop()

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
        self.combobox = ttkbs.Combobox(
            master=self.Frame1,
            bootstyle=DANGER,
            font=('微软雅黑',10),
            values=['blastn','tblastx'],
            width=10,textvariable=self.combobox_var
        )
        self.combobox.current(0)
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
        self.combobox.grid(row=3,column=0)
        self.Entry_user.grid(row=1, column=2, pady=3)
        self.Entry_password.grid(row=2, column=2, pady=4)
        self.Entry_parameters.grid(row=3,column=2,pady=4)
        self.button_go.grid(row=4, column=2, pady=9)

    # cmd指令调整
    def blast(self,way,db_file,blast_file,c_db,parameters=''):
        self.ret = {}
        cmd1 = f'makeblastdb -in {db_file} -dbtype nucl -out {db_file}.blastdb'
        cmd2 = f'{way} -query {blast_file} -db {db_file}.blastdb -out {c_db} -outfmt 6 {parameters}'
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
            name = new_name + str(self.num-2) + fa
            self.Flag2 = os.path.getsize(f'./{name}')

    # 修改文件操作
    def clean(self,file_c,file_a):
        # 创建新的B库
        with open(file_c, 'w', encoding='utf-8') as fp:
            for i in self.name:
                if i in self.data_dict.keys():
                    fp.write(self.data_dict[i])
                    self.save_data[i] = self.data_dict[i]
                    del self.data_dict[i]

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

    # 制作数据库名称函数
    def db_name_maker(self):
        self.c_db = self.db_file
        c_db, fa = os.path.splitext(self.c_db)
        self.c_db = c_db + str(self.num) + '.blast'

    #给blast数据库加数字
    def c_blast_name(self):
        c_db, fa = os.path.splitext(self.c_db)
        self.c_db = c_db + str(self.num) + fa

    # 替换R|Y|M|K|S|W|H|B|V|D特殊字符为N的函数
    def sub(self,file_name):
        with open(f"{file_name}", encoding='utf-8') as file:
            strings = file.read()
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
        self.way = self.combobox_var.get()

        # 开始流程
        try:
            # 替换所有R|Y|M|K|S|W|H|B|V|D
            self.sub(self.blast_file)
            self.sub(self.db_file)
            while True:
                # 判断有无参数，并执行cmd语句
                if self.parameters:
                    self.blast(self.way, self.db_file, self.blast_file, self.c_db,parameters=self.parameters)
                else:
                    self.blast(self.way, self.db_file, self.blast_file,self.c_db)
                # 先读数据（便于后续对比）
                self.read_file()
                if self.Flag1 == 0 or self.Flag2 == self.Flag1:
                    break
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
            self.label_var.set('Blast Successfully')
            messagebox.showinfo('Success','Blast Successfully')
            # 创建A_iterations文件
            with open('A_iterations.txt','w',encoding='utf-8') as fp:
                fp.write(f'count:{self.count-1}')
        except Exception as e:
            messagebox.showinfo('Error',e)
            self.label_var.set('')
if __name__ == '__main__':
    Blast = Blast()