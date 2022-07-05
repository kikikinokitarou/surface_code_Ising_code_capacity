import openjij as oj
from pyqubo import Array, Constraint, Placeholder, solve_qubo , Xor,Binary
import numpy as np
from numpy import zeros
import re
import copy
import random
import cxxjij.utility as U
import math
import time






ls=7 #lattice size
num_qubit=ls*ls+(ls-1)*(ls-1)       #number of qubit
num_stabilizer=ls*(ls-1)        #number of stabilizer
syn =[ [0 for j in range(num_stabilizer)] for i in range(2)]        #シンドローム
stabilizer=[[[0 for i in range(4)] for j in range(num_stabilizer)]for k in range(2)]        #スタビライザー用の行列(後程生成)
qubit=[ [0 for j in range(num_qubit)] for i in range(2)]        #量子ビットのエラー情報
result=[[0 for j in range(num_stabilizer)] for i in range(2)]       #SAによる解
ER=[ [0 for j in range(num_qubit)] for i in range(3)]       #推定したエラー
ER_x=[ [0 for j in range(num_qubit)] for i in range(3)]       #推定したエラー
ER_z=[ [0 for j in range(num_qubit)] for i in range(3)]       #推定したエラー
ER_xz=[ [0 for j in range(num_qubit)] for i in range(3)]       #推定したエラー
logical=[[0 for j in range(ls)] for i in range(2)]       #logical error判定用の行列
start=[0 for i in range(2*num_stabilizer)]

list_p_logi=[]
list_errorbar=[]
list_p_phys=[]
list_excution_time=[]


#バイナリー変数
x=Array.create('x',shape=(num_stabilizer), vartype='BINARY')
z=Array.create('z',shape=(num_stabilizer), vartype='BINARY')

num_sample=10000
    #number of sample




#(0,1)から(1,-1)へ
def binary_to_spin(b):
    return 1-2*b

#(1,-1)から(1,0)へ
def spin_to_binary(s):
    return (s+1)/2

def gene_logical():
    for i in range(ls):
        logical[0][i]=i
        logical[1][i]=ls*i


#スタビライザーの生成
def gene_sta():
    for i in range(num_stabilizer):
        stabilizer[0][i][0]=i
        stabilizer[0][i][1]=i+ls
        if i%ls==0:
            stabilizer[0][i][2]=i + ls**2 - i//ls
            stabilizer[0][i][3]=stabilizer[0][i][2]
        elif i%ls==ls-1:
            stabilizer[0][i][2]=i + ls**2 - i//ls - 1
            stabilizer[0][i][3]=stabilizer[0][i][2]
        else:
            stabilizer[0][i][2]=i + ls**2 - i//ls - 1
            stabilizer[0][i][3]=stabilizer[0][i][2]+1
    for i in range(num_stabilizer):
        stabilizer[1][i][0]=i + i//(ls-1)
        stabilizer[1][i][1]=stabilizer[1][i][0] + 1
        if 0<=i<=ls-2:
            stabilizer[1][i][2]=i + ls**2
            stabilizer[1][i][3]=stabilizer[1][i][2]
        elif (ls-1)**2<=i<=ls*(ls-1)-1:
            stabilizer[1][i][2]=i + ls**2 - ls + 1
            stabilizer[1][i][3]=stabilizer[1][i][2]
        else:
            stabilizer[1][i][3]=i + ls**2
            stabilizer[1][i][2]=stabilizer[1][i][3] - ls + 1

def gene_SAschedule(p):
    beta=-math.log(p/(1-p))/2
    schedule_list = U.make_classical_schedule_list(0.0001, beta, num_qubit, 100)
    return schedule_list


def initialize():
    for i in range(num_qubit):
        qubit[0][i]=0
        qubit[1][i]=0


def gene_error(p_x,p_z,p_y):
    for i in range(num_qubit):
        p=random.random()
        if p<p_x:
            qubit[0][i]^=1
        elif p<(p_x+p_z):
            qubit[1][i]^=1
        elif p<(p_x+p_z+p_y):
            qubit[0][i]^=1
            qubit[1][i]^=1







#シンドローム測定
def syndrome_meas():
    for i in range(2):
        for j in range(num_stabilizer):
            syn[i][j]=0
            for k in range(4):
                if k==0:
                    syn[i][j]^=qubit[i][stabilizer[i][j][k]]
                elif stabilizer[i][j][k]!=stabilizer[i][j][k-1]:
                    syn[i][j]^=qubit[i][stabilizer[i][j][k]]

#シンドローム測定結果の確認
def show_syndrome():
    for i in range(2):
        for j in range(num_stabilizer):
            if(syn[i][j]==1):
                if(i==0):
                    print('x:',j)
                if(i==1):
                    print("z: ",j)

#初期エラーの生成
def gene_pureerror():
    ER_pure=[ [0 for j in range(num_qubit)] for i in range(2)]      #エラーの初期配置

    #Xエラーについて
    for i in range(num_stabilizer):
        if(i<num_stabilizer//2):
            if syn[0][i]==1:
                for j in range(i//ls+1):
                    ER_pure[0][i-ls*j]^=1
        else:
            if syn[0][i]==1:
                for j in range(1,ls-i//ls):
                    ER_pure[0][i+ls*j]^=1

    #Zエラーについて
    for i in range(num_stabilizer):
        if(i%(ls-1)<(ls-1)//2):
            if syn[1][i]==1:
                for j in range(i%(ls-1)+1):
                    ER_pure[1][ls*(i//(ls-1))+j]^=1
        else:
            if syn[1][i]==1:
                for j in range((ls-1)-i%(ls-1)):
                    ER_pure[1][i+i//(ls-1)+1+j]^=1
    return ER_pure
    

def gene_pureerror_x(ER_pure):
    #初期エラーの更新
    ER_pure_x=copy.deepcopy(ER_pure)
    for i in range(ls):
        ER_pure_x[0][ls*i]^=1
    return ER_pure_x

def gene_pureerror_z(ER_pure):
    #初期エラーの更新
    ER_pure_z=copy.deepcopy(ER_pure)
    for i in range(ls):
        ER_pure_z[1][i]^=1
    return ER_pure_z

def gene_pureerror_xz(ER_pure):
    #初期エラーの更新
    ER_pure_xz=copy.deepcopy(ER_pure)
    for i in range(ls):
        ER_pure_xz[1][i]^=1
        ER_pure_xz[0][ls*i]^=1
    return ER_pure_xz

#コスト関数の作り方
def gene_costfunc(ER_pure):

    H_cost=0
    #Xエラーについて
    for i in range(num_qubit):
        if i<ls**2:
            if i%ls==0:
                H_cost-=1*binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls])
            elif i%ls==ls-1:
                H_cost-=1*binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls-1])
            else:
                H_cost-=1*binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls-1])*binary_to_spin(x[i-i//ls])
        else:
            H_cost-=1*binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-ls**2])*binary_to_spin(x[i-ls**2+ls-1])
            

    #Zエラーについて
    for i in range(num_qubit):
        if i<ls**2:
            if i<ls:
                H_cost-=1*binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i])
            elif i>=ls*(ls-1):
                H_cost-=1*binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i-ls])
            else:
                H_cost-=1*binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i-ls])*binary_to_spin(z[i])
        else:
            H_cost-=1*binary_to_spin(ER_pure[0][i])*binary_to_spin(z[i-ls**2+(i-ls**2)//(ls-1)])*binary_to_spin(z[i-ls**2+(i-ls**2)//(ls-1)+1])
            
        
    #Yエラーについて(XエラーとZエラーを打ち消す)

    #2つのスタビライザーに挟まれたビット
    H_cost-=2*spin_to_binary( -binary_to_spin(ER_pure[0][0])*binary_to_spin(x[0]) )*spin_to_binary( -binary_to_spin(ER_pure[1][0])*binary_to_spin(z[0]) )
    H_cost-=2*spin_to_binary(- binary_to_spin(ER_pure[0][ls-1])*binary_to_spin(x[ls-2]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls-1])*binary_to_spin(z[ls-1]) )
    H_cost-=2*spin_to_binary(-binary_to_spin(ER_pure[0][ls*(ls-1)])*binary_to_spin(x[(ls-1)**2]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*(ls-1)])*binary_to_spin(z[ls*(ls-2)]))
    H_cost-=2*spin_to_binary(- binary_to_spin(ER_pure[0][ls**2-1])*binary_to_spin(x[ls*(ls-1)-1]))*spin_to_binary(- binary_to_spin(ER_pure[1][ls**2-1])*binary_to_spin(z[ls*(ls-1)-1]))

    #3つのスタビライザーに囲まれた量子ビット
    for i in range(1,ls-1):
        H_cost-=2*spin_to_binary(-binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls-1])*binary_to_spin(x[i-i//ls]))*spin_to_binary(-binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i]))
        H_cost-=2*spin_to_binary(-binary_to_spin(ER_pure[0][ls*i])*binary_to_spin(x[(ls-1)*i]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*i])*binary_to_spin(z[ls*i-ls])*binary_to_spin(z[ls*i]))
        H_cost-=2*spin_to_binary(-binary_to_spin(ER_pure[0][ls*i+(ls-1)])*binary_to_spin(x[ls*i+(ls-1)-(i+1)]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*i+(ls-1)])*binary_to_spin(z[ls*i+(ls-1)-ls])*binary_to_spin(z[ls*i+(ls-1)]))
        H_cost-=2*spin_to_binary(-binary_to_spin(ER_pure[0][ls*(ls-1)+i])*binary_to_spin(x[ls*(ls-1)+i-ls])*binary_to_spin(x[ls*(ls-1)+i-ls+1]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*(ls-1)+i])*binary_to_spin(z[ls*(ls-1)+i-ls]))
        
        
    #4つのスタビライザーに囲まれた量子ビット
    for i in range(1,ls-1):
        for j in range(1,ls-1):
            H_cost-=2*spin_to_binary(-binary_to_spin(ER_pure[0][ls*i+j])*binary_to_spin(x[(ls-1)*i+j-1])*binary_to_spin(x[(ls-1)*i+j]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*i+j])*binary_to_spin(z[ls*(i-1)+j])*binary_to_spin(z[ls*i+j]))

            
            
    for i in range(ls-1):
        for j in range(ls-1):
            H_cost-=2*spin_to_binary(-binary_to_spin(ER_pure[0][ls**2+(ls-1)*i+j])*binary_to_spin(x[(ls-1)*i+j])*binary_to_spin(x[(ls-1)*(i+1)+j]))*spin_to_binary(-binary_to_spin(ER_pure[0][ls**2+(ls-1)*i+j])*binary_to_spin(z[ls*i+j])*binary_to_spin(z[ls*i+j+1]))


    #Yエラーについて（Yエラーを数える）

    #2つのスタビライザーに挟まれたビット
    H_cost+=1*spin_to_binary( -binary_to_spin(ER_pure[0][0])*binary_to_spin(x[0]) )*spin_to_binary( -binary_to_spin(ER_pure[1][0])*binary_to_spin(z[0]) )
    H_cost+=1*spin_to_binary(- binary_to_spin(ER_pure[0][ls-1])*binary_to_spin(x[ls-2]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls-1])*binary_to_spin(z[ls-1]) )
    H_cost+=1*spin_to_binary(-binary_to_spin(ER_pure[0][ls*(ls-1)])*binary_to_spin(x[(ls-1)**2]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*(ls-1)])*binary_to_spin(z[ls*(ls-2)]))
    H_cost+=1*spin_to_binary(- binary_to_spin(ER_pure[0][ls**2-1])*binary_to_spin(x[ls*(ls-1)-1]))*spin_to_binary(- binary_to_spin(ER_pure[1][ls**2-1])*binary_to_spin(z[ls*(ls-1)-1]))

    #3つのスタビライザーに囲まれた量子ビット
    for i in range(1,ls-1):
        H_cost+=1*spin_to_binary(-binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls-1])*binary_to_spin(x[i-i//ls]))*spin_to_binary(-binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i]))
        H_cost+=1*spin_to_binary(-binary_to_spin(ER_pure[0][ls*i])*binary_to_spin(x[(ls-1)*i]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*i])*binary_to_spin(z[ls*i-ls])*binary_to_spin(z[ls*i]))
        H_cost+=1*spin_to_binary(-binary_to_spin(ER_pure[0][ls*i+(ls-1)])*binary_to_spin(x[ls*i+(ls-1)-(i+1)]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*i+(ls-1)])*binary_to_spin(z[ls*i+(ls-1)-ls])*binary_to_spin(z[ls*i+(ls-1)]))
        H_cost+=1*spin_to_binary(-binary_to_spin(ER_pure[0][ls*(ls-1)+i])*binary_to_spin(x[ls*(ls-1)+i-ls])*binary_to_spin(x[ls*(ls-1)+i-ls+1]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*(ls-1)+i])*binary_to_spin(z[ls*(ls-1)+i-ls]))
        
        
    #4つのスタビライザーに囲まれた量子ビット
    for i in range(1,ls-1):
        for j in range(1,ls-1):
            H_cost+=1*spin_to_binary(-binary_to_spin(ER_pure[0][ls*i+j])*binary_to_spin(x[(ls-1)*i+j-1])*binary_to_spin(x[(ls-1)*i+j]))*spin_to_binary(-binary_to_spin(ER_pure[1][ls*i+j])*binary_to_spin(z[ls*(i-1)+j])*binary_to_spin(z[ls*i+j]))

            
            
    for i in range(ls-1):
        for j in range(ls-1):
            H_cost+=1*spin_to_binary(-binary_to_spin(ER_pure[0][ls**2+(ls-1)*i+j])*binary_to_spin(x[(ls-1)*i+j])*binary_to_spin(x[(ls-1)*(i+1)+j]))*spin_to_binary(-binary_to_spin(ER_pure[0][ls**2+(ls-1)*i+j])*binary_to_spin(z[ls*i+j])*binary_to_spin(z[ls*i+j+1]))
    return H_cost*10

#最適化
def solve(H,schedule_list):
    model=H.compile(strength=30)
    qubo , offset=model.to_qubo()
    sampler=oj.SASampler(start,schedule=schedule_list)
    response=sampler.sample_qubo(qubo)
    for i in range(2):
        for j in range(num_stabilizer):
            result[i][j]=0
    for i in range(2):
        for j in range(num_stabilizer):
            if i==0:
                s="x["+str(j)+"]"
                result[i][j]=response.first.sample[s]
            else:
                s='z['+str(j)+']'
                result[i][j]=response.first.sample[s]



#エラーの推定結果
def est_error(ER_pure):
    error=[ [0 for j in range(num_qubit)] for i in range(3)]

    for i in range(num_qubit):
        if i<ls**2:
            if i%ls==0:
                error[0][i]^=ER_pure[0][i]^result[0][i-i//ls]
            elif i%ls==ls-1:
                error[0][i]^=ER_pure[0][i]^result[0][i-i//ls-1]
            else:
                error[0][i]^=ER_pure[0][i]^result[0][i-i//ls-1]^result[0][i-i//ls]
        else:
            error[0][i]^=ER_pure[0][i]^result[0][i-ls**2]^result[0][i-ls**2+ls-1]
            
    for i in range(num_qubit):
        if i<ls**2:
            if i<ls:
                error[1][i]^=ER_pure[1][i]^result[1][i]
            elif i>=ls*(ls-1):
                error[1][i]^=ER_pure[1][i]^result[1][i-ls]
            else:
                error[1][i]^=ER_pure[1][i]^result[1][i-ls]^result[1][i]
        else:
            error[1][i]^=ER_pure[0][i]^result[1][i-ls**2+(i-ls**2)//(ls-1)]^result[1][i-ls**2+(i-ls**2)//(ls-1)+1]
        
            
    for i in range(num_qubit):
        if(error[0][i]==1 and error[1][i]==1):
            error[2][i]^=1
            error[0][i]^=1
            error[1][i]^=1
    return error
    
    

#推定したエラーの出力
def show_result(error):
    for i in range(num_qubit):
        if(error[0][i]==1):
            print('x: ',i)
    for i in range(num_qubit):
        if(error[1][i]==1):
            print('z: ',i)
    for i in range(num_qubit):
        if(error[2][i]==1):
            print('y: ',i)

#推定したエラーの出力
def show_error(error):
    for i in range(num_qubit):
        if(error[0][i]==1):
            print('x: ',i)
    for i in range(num_qubit):
        if(error[1][i]==1):
            print('z: ',i)
#発生したエラーの個数のカウント

#推定したエラーの個数のカウント
def count_result(error):
    l=0
    for i in range(num_qubit):
            l+=error[0][i]+error[1][i]
    return l

def recovery(error):
    for i in range(num_qubit):
        qubit[0][i]^=error[0][i]^error[2][i]
        qubit[1][i]^=error[1][i]^error[2][i]

def judge_logical():
    flg_logi_x=0
    flg_logi_z=0
    for i in range(ls):
        flg_logi_x ^= qubit[0][logical[0][i]]
        flg_logi_z^=qubit[1][logical[1][i]]
    return (flg_logi_x!=0 or flg_logi_z!=0)
                






gene_sta()
gene_logical()


for p in range(3,11):
    p_phys=p/100
    p_phys+=0.005
    p_logi=0
    events_logi=0
    schedule_list=gene_SAschedule(2*p_phys/3)
    excution_time=0
    for i in range(3):
        initialize()
        gene_error(p_phys/3,p_phys/3,p_phys/3)

        syndrome_meas()

        start1=time.time()
        ER_pure=gene_pureerror()
        end1=time.time()


        H=gene_costfunc(ER_pure)
        solve(H,schedule_list)
        ER=est_error(ER_pure)
        l=count_result(ER)

        start2=time.time()
        ER_pure_x=gene_pureerror_x(ER_pure)
        H_x=gene_costfunc(ER_pure_x)
        solve(H_x,schedule_list)
        ER_x=est_error(ER_pure_x)
        l_x=count_result(ER_x)
        end2=time.time()


        ER_pure_z=gene_pureerror_z(ER_pure)
        H_z=gene_costfunc(ER_pure_z)
        solve(H_z,schedule_list)
        ER_z=est_error(ER_pure_z)
        l_z=count_result(ER_z)


        ER_pure_xz=gene_pureerror_xz(ER_pure)
        H_xz=gene_costfunc(ER_pure_xz)
        solve(H_xz,schedule_list)
        ER_xz=est_error(ER_pure_xz)
        l_xz=count_result(ER_xz)


        l_min=min(l,l_x,l_z,l_xz)
    
        excution_time+=end1-start1+end2-start2
        #print("発生したエラー")
        #Sshow_error(qubit)
        if(l_min==l):
            #print("初期エラー")
            #show_error(ER_pure)
            #print("推定したエラー(論理エラーなし)")
            #show_result(ER)
            #print("初期エラー(論理エラーX)")
            #show_error(ER_pure_x)
            #print("推定したエラー(論理エラーX)")
            #show_result(ER_x)
            #print("初期エラー(論理エラーZ)")
            #show_error(ER_pure_z)
            #print("推定したエラー(論理エラーZ)")
            #show_result(ER_z)
            recovery(ER)
            #print("訂正後")
            #show_error(qubit)
            #print("エラーチェーンの長さ")
            #print(l,l_x,l_z,l_xz)
        elif(l_min==l_x):
            #print("初期エラー(論理エラーなし)")
            #show_error(ER_pure)
            #print("初期エラー")
            #show_error(ER_pure_x)
            #print("推定したエラー(論理エラーX)")
            #show_result(ER_x)
            recovery(ER_x)
            #print("訂正後")
            #show_error(qubit)
            #print("エラーチェーンの長さ")
            #print(l,l_x,l_z,l_xz)
        elif(l_min==l_z):
            #print("初期エラー(論理エラーなし)")
            #show_error(ER_pure)
            #print("初期エラー")
            #show_error(ER_pure_z)
            #print("推定したエラー(論理エラーZ)")
            #show_result(ER_z)
            recovery(ER_z)
            #print("訂正後")
            #show_error(qubit)
            #print("エラーチェーンの長さ")
            #print(l,l_x,l_z,l_xz)
        elif(l_min==l_xz):
            #print("初期エラー(論理エラーなし)")
            #show_error(ER_pure)
            #print("初期エラー")
            #show_error(ER_pure_xz)
            #print("推定したエラー(論理エラーXZ)")
            #show_error(ER_xz)
            recovery(ER_xz)
            #print("訂正後")
            #show_error(qubit)
            #print("エラーチェーンの長さ")
            #print(l,l_x,l_z,l_xz)
        
        if(judge_logical()):
            print("logical_error")
            #print(p,i)
            #exit()
            events_logi+=1
        print(p_phys,i)
    p_logi=events_logi/num_sample
    errorbar=math.sqrt(events_logi)/num_sample
    list_p_phys.append(p_phys)
    list_p_logi.append(p_logi)
    list_errorbar.append(errorbar)
    list_excution_time.append(excution_time/num_sample)

print(list_p_logi)
print(list_errorbar)
print(list_excution_time)





        


        
    











                    

