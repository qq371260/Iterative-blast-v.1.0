from subprocess import Popen,PIPE
import os
import re
import tkinter as tk
from ttkbootstrap.constants import *
import ttkbootstrap as ttkbs
from tkinter import messagebox
class Blast(object):
    def __init__(self):
        self.root = ttkbs.Window()
        self.root.title('Blast')
        self.root.geometry('550x400+400+200')

        # Variable reception
        self.user_var = ttkbs.StringVar()
        self.password_var = ttkbs.StringVar()
        self.parameters_var = ttkbs.StringVar()
        self.label_var = ttkbs.StringVar()
        self.user_var.set('B.fasta.fa')
        self.password_var.set('A.fasta.fa')
        self.combobox_var = ttkbs.StringVar()
        self.combobox_var.set('')
        self.num = 1

        # Data reception
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

    def setting(self):

        # Page setup
        self.Frame1 = tk.Frame(self.root, background='black', bg='blue')

        # GUISettings
        self.blank = ttkbs.Label(self.Frame1, width=4, border=10)
        self.Label1 = ttkbs.Label(self.Frame1, text='db file:')
        self.Label2 = ttkbs.Label(self.Frame1, text='query file:')
        self.Label3 = ttkbs.Label(self.Frame1, text='parameters:')
        self.Label_notice = ttkbs.Label(self.Frame1,font=('Msyh',12),textvariable=self.label_var)
        self.combobox = ttkbs.Combobox(
            master=self.Frame1,
            bootstyle=DANGER,
            font=('Msyh',10),
            values=['blastn','blastx','tblastn','blastp'],
            width=10,textvariable=self.combobox_var
        )
        self.combobox.current(0)
        self.Entry_user = ttkbs.Entry(self.Frame1, textvariable=self.user_var)
        self.Entry_password = ttkbs.Entry(self.Frame1, textvariable=self.password_var)
        self.Entry_parameters = ttkbs.Entry(self.Frame1,textvariable=self.parameters_var)

        self.button_go = ttkbs.Button(self.Frame1, text='Go', bootstyle=PRIMARY, width=10,command=self.go)

        # Grid layouts
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
    def blast(self,way,db_file,blast_file,c_db,parameters=''):
        self.ret = {}
        cmd1 = f'makeblastdb -in {db_file} -dbtype nucl -out {db_file}.blastdb'
        cmd2 = f'{way} -query {blast_file} -db {db_file}.blastdb -out {c_db} -outfmt 6 -evalue 1e-{parameters}'
        p = Popen(cmd1 + '&&' + cmd2, stdin=PIPE, stdout=PIPE,
                             stderr=PIPE, shell=True)
        p.wait()
    def read_file(self):
        # Operating B library
        new_name = self.c_db + '.txt'
        os.rename(self.c_db, new_name)
        self.c_db = new_name
        with open(self.c_db, 'r', encoding='utf-8') as fp:
            content = fp.read()
            self.Flag1 = content
            columns = re.findall(r'(^\w.*?)\t', content, re.MULTILINE)
            # print(columns)
            for column in columns:
                self.name.add(column)
        c_db, txt = os.path.splitext(self.c_db)
        os.rename(self.c_db,c_db)
        self.c_db = c_db

        # Operating A library
        self.Flag1 = os.path.getsize(f'./{self.blast_file}')
        self.rename_txt()
        with open(self.blast_file, 'r', encoding='utf-8') as fp:
            data = fp.read()
            name_pre = re.findall(r'^>(.*)', data, re.MULTILINE)
            result = re.findall(r'>.[^>]*', data, re.S)
            self.data_dict = {name_pre: result for name_pre, result in zip(name_pre, result)}

        self.rename_fa()
        name_fasta, fa = os.path.splitext(str(self.blast_file))
        new_name, fasta = os.path.splitext(name_fasta)
        int_ = re.compile('\d')
        new_name = re.sub(int_,'', new_name)
        if self.num >= 3:
            name = new_name + str(self.num-2) + fasta + fa
            self.Flag2 = os.path.getsize(f'./{name}')

    def clean(self,file_c,file_a):
        # Creating new B library
        with open(file_c, 'w', encoding='utf-8') as fp:
            for i in self.name:
                if i in self.data_dict.keys():
                    fp.write(self.data_dict[i])
                    self.save_data[i] = self.data_dict[i]
                    del self.data_dict[i]

        # Creating new A library
        with open(file_a, 'w', encoding='utf-8') as fp:
            for i in self.data_dict.values():
                fp.write(i)
    #Only A file
    def rename_txt(self):
        new_name = self.blast_file + '.txt'
        os.rename(self.blast_file,new_name)
        self.blast_file = new_name

    def rename_fa(self):
        new_name,txt = os.path.splitext(self.blast_file)
        os.rename(self.blast_file,new_name)
        self.blast_file = new_name

    def add_file_num(self,name):
        name_fasta, fa = os.path.splitext(name)
        new_name,fasta = os.path.splitext(name_fasta)
        int_ = re.compile('\d')
        new_name = re.sub(int_,"",new_name)
        name = new_name + str(self.num) + fasta + fa
        return name
    def db_name_maker(self):
        self.c_db = self.db_file
        c_db_fasta, fa = os.path.splitext(self.c_db)
        c_db, fasta = os.path.splitext(c_db_fasta)
        self.c_db = c_db + str(self.num) + '.blast'
    def c_blast_name(self):
        c_db_fasta, fa = os.path.splitext(self.c_db)
        c_db, fasta = os.path.splitext(c_db_fasta)
        self.c_db = c_db + str(self.num) + fasta + fa
    def go(self):
        self.db_file =self.user_var.get()
        self.blast_file = self.password_var.get()
        self.db_name_maker()
        self.parameters = self.parameters_var.get()
        self.way = self.combobox_var.get()

        # Beginning
        try:
            while  self.Flag1 != 0 and self.Flag2 != self.Flag1:
                # Checking parameters and excuting CMD statement
                if self.parameters:
                    self.blast(self.way, self.db_file, self.blast_file, self.c_db,parameters=self.parameters)
                    # print('sb')
                else:
                    self.blast(self.way, self.db_file, self.blast_file,self.c_db)
                    # print('cg')
                # Reading data for later BLAST
                self.read_file()
                # Adding number
                self.db_file = self.add_file_num(self.db_file)
                self.blast_file = self.add_file_num(self.blast_file)
                # Data clean, genetating new .fa file
                self.clean(self.db_file,self.blast_file)
                self.num += 1
                self.count += 1
            # Execution end
            with open('B_total.txt','w',encoding='utf-8') as fp:
                for data in self.save_data.values():
                    fp.write(data)
            self.label_var.set('Blast Successfully')
            messagebox.showinfo('Success','Blast Successfully')
            with open('A_iterations.txt','w',encoding='utf-8') as fp:
                fp.write(f'count:{self.count-2}')
                a = 'A_iterations.txt'
                name , txt = os.path.splitext(a)
                new = name + '.fasta.fa'
            os.rename(a,new)
        except Exception as e:
            messagebox.showinfo('error','An error occurred')
            self.label_var.set('')
            print(e)
            self.label_var.set('Blast Successfully')
if __name__ == '__main__':
    Blast = Blast()
