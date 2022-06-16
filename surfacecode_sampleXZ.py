import openjij as oj
from pyqubo import Array, Constraint, Placeholder, solve_qubo , Xor,Binary
import numpy as np
from numpy import zeros
import re
import copy
import sys
import pprint

""""

注意：このサンプルファイルには初期エラーに論理演算子を含むパターンを導入していない


"""


ls=5 #lattice size
num_qubit=ls*ls+(ls-1)*(ls-1)       #number of qubit
num_stabilizer=ls*(ls-1)        #number of stabilizer
syn =[ [0 for j in range(num_stabilizer)] for i in range(2)]        #シンドローム
stabilizer=[[[0 for i in range(4)] for j in range(num_stabilizer)]for k in range(2)]        #スタビライザー用の行列(後程生成)
qubit=[ [0 for j in range(num_qubit)] for i in range(2)]        #量子ビットのエラー情報
ER_pure=[ [0 for j in range(num_qubit)] for i in range(2)]      #エラーの初期配置
result=[[0 for j in range(num_stabilizer)] for i in range(2)]       #SAによる解
ER=[ [0 for j in range(num_qubit)] for i in range(3)]       #推定したエラー

#バイナリー変数
x=Array.create('x',shape=(num_stabilizer), vartype='BINARY')
z=Array.create('z',shape=(num_stabilizer), vartype='BINARY')



#(0,1)から(1,-1)へ
def binary_to_spin(b):
    return 1-2*b

#(1,-1)から(1,0)へ
def spin_to_binary(s):
    return (s+1)/2


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
def gene_pureerror(ER_pure):
    for i in range(2):
        for j in range(num_qubit):
            ER_pure[i][j]=0

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

#コスト関数の作り方
def gene_costfunc():

    H_cost=0
    #Xエラーについて
    for i in range(num_qubit):
        if i<ls**2:
            if i%ls==0:
                H_cost-=binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls])
            elif i%ls==ls-1:
                H_cost-=binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls-1])
            else:
                H_cost-=binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-i//ls-1])*binary_to_spin(x[i-i//ls])
        else:
            H_cost-=binary_to_spin(ER_pure[0][i])*binary_to_spin(x[i-ls**2])*binary_to_spin(x[i-ls**2+ls-1])
            

    #Zエラーについて
    for i in range(num_qubit):
        if i<ls**2:
            if i<ls:
                H_cost-=binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i])
            elif i>=ls*(ls-1):
                H_cost-=binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i-ls])
            else:
                H_cost-=binary_to_spin(ER_pure[1][i])*binary_to_spin(z[i-ls])*binary_to_spin(z[i])
        else:
            H_cost-=binary_to_spin(ER_pure[0][i])*binary_to_spin(z[i-ls**2+(i-ls**2)//(ls-1)])*binary_to_spin(z[i-ls**2+(i-ls**2)//(ls-1)+1])
            
        
    return H_cost

#最適化
def solve(H):
    model=H.compile()
    qubo , offset=model.to_qubo()
    sampler=oj.SASampler(num_reads=1,num_sweeps=100000)
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
def est_error(ER):
    for i in range(3):
        for j in range(num_qubit):
            ER[i][j]=0

    for i in range(num_qubit):
        if i<ls**2:
            if i%ls==0:
                ER[0][i]^=ER_pure[0][i]^result[0][i-i//ls]
            elif i%ls==ls-1:
                ER[0][i]^=ER_pure[0][i]^result[0][i-i//ls-1]
            else:
                ER[0][i]^=ER_pure[0][i]^result[0][i-i//ls-1]^result[0][i-i//ls]
        else:
            ER[0][i]^=ER_pure[0][i]^result[0][i-ls**2]^result[0][i-ls**2+ls-1]
            
    for i in range(num_qubit):
        if i<ls**2:
            if i<ls:
                ER[1][i]^=ER_pure[1][i]^result[1][i]
            elif i>=ls*(ls-1):
                ER[1][i]^=ER_pure[1][i]^result[1][i-ls]
            else:
                ER[1][i]^=ER_pure[1][i]^result[1][i-ls]^result[1][i]
        else:
            ER[1][i]^=ER_pure[0][i]^result[1][i-ls**2+(i-ls**2)//(ls-1)]^result[1][i-ls**2+(i-ls**2)//(ls-1)+1]
        
            
    for i in range(num_qubit):
        if(ER[0][i]==1 and ER[1][i]==1):
            ER[2][i]^=1
            ER[0][i]^=1
            ER[1][i]^=1
    return ER
    

#推定したエラーの出力
def show_error(ER):
    for i in range(num_qubit):
        if(ER[0][i]==1):
            print('x: ',i)
    for i in range(num_qubit):
        if(ER[1][i]==1):
            print('z: ',i)
    for i in range(num_qubit):
        if(ER[2][i]==1):
            print('y: ',i)
#エラーの個数のカウント
def count_error(ER):
    l=0
    for i in range(num_qubit):
        if(ER[0][i]==1 or ER[1][i]==1 or ER[2][i]==1):
            l+=1
    return l


#実際に生じるエラー
qubit[0][5]=1
qubit[0][29]=1
qubit[1][1]=1
qubit[1][26]=1

gene_sta()
syndrome_meas()
#show_syndrome()
ER_pure=gene_pureerror(ER_pure)
cou=0
for i in range(1000):
    H=gene_costfunc()
    solve(H)
    est_error(ER)
    l=count_error(ER)
    if(l>2):
        cou+=1
    print(i)
print(cou)
    











                    

