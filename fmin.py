#!/use/bin/env python
# -*- coding: UTF-8 -*-

from numpy import *
from scipy import *
from scipy import linalg
import string
import sys

def main():
    args = sys.argv[1:]
    problem = int(args[0]) #課題を決める変数
    method = int(args[1]) #手法を決める変数
    search = int(args[2]) #直線探索を決める変数

    #各課題についてevalf,evalg,evalhを定義,初期点(x)・次元数(n)を決定
    if problem == 1:
        n = 2
        x = zeros(n)
        for i in range(0,n): #初期点を入力
            x[i] = float(input('初期点の第'+str(i+1)+'成分を入力してください: '))
        x_0 = x

        def evalf(x):
            f = x[0]**2 + math.exp(x[0]) + x[1]**4 + x[1]**2 - 2*x[0]*x[1] + 3
            return f
        def evalg(x):
            g = array ([2*x[0] +  math.exp(x[0]) - 2*x[1],4*x[1]**3 + 2*x[1] - 2*x[0]])
            return g
        def evalh(x):
            H = array([[2 + math.exp(x[0]),-2],[-2,12*x[1]**2 + 2]])
            return H

    elif problem == 2:
        n = 2
        x = array([1.,1.])

        def evalf(x):
            f = 0.0
            for i in range(0,3):
                f += evalf_i(i,x)**2
            return f

        def evalf_i(i,x):
            y = array([1.5,2.25,2.625])
            f = y[i] - x[0] * (1 - x[1]**(i+1))
            return f
        def evalg_i(i,x):
            g = array([x[1]**(i+1) - 1,(i+1) * x[0] * x[1]**i])
            return g
        def evalh_i(i,x):
            h = array([[0,(i+1) * x[1]**i],[(i+1) * x[1]**i,(i+1) * i * x[0] * x[1]**(i-1)]])
            return h

        def evalg(x):
            g = zeros(2)
            for i in range(0,3):
                g += 2 * evalf_i(i,x) * evalg_i(i,x)
            return g
        def evalh(x):
            H = array([[0.0,0.0],[0.0,0.0]])
            for i in range(0,3):
                H += 2 * (evalf_i(i,x) * evalh_i(i,x) + dot(evalg_i(i,x).reshape(-1,1),evalg_i(i,x).reshape(1,-1)))
            return H

    elif problem == 3:
        n = int(input('次元数: ')) #次元数を入力

        def evalf(x):
            f = inner(dot(A,x),x) / 2.0
            return f
        def evalg(x):
            g = dot(A,x)
            return g
        def evalh(x):
            H = A
            return H

    elif problem == 4:
        n = 2
        x = array([0.,0.])
        l = array([-5,0]) #制約集合の下限
        u = array([10,15]) #制約集合の上限

        def evalf(x):
            f = (x[1] - 5.1*x[0]**2 / (4*math.pi**2) + 5*x[0] / math.pi -6)**2 + 10*(1 - 1/(8*math.pi))*math.cos(x[0]) + 10
            return f
        def evalg(x):
            g = array([2*(x[1] - 5.1*x[0]**2 / (4*math.pi**2) + 5*x[0] / math.pi -6)*(-5.1*x[0] / (2*math.pi**2) + 5/math.pi) - 10*(1 - 1/ (8*math.pi))*math.sin(x[0]),2*(x[1] - 5.1*x[0]**2 / (4*math.pi**2) + 5*x[0] / math.pi -6)])
            return g

    elif problem == 5:
        n = int(input('次元数: ')) #次元数を入力
        x = zeros(n)
        for i in range(0,n): #初期点を入力
          x[i] = float(input('初期点の第'+str(i+1)+'成分を入力してください: '))
        x_0 = x

        l = -32.768 * ones(n) #制約集合の下限
        u = 32.768 * ones(n) #制約集合の上限

        def evalf1(x):
            f1 = 0
            for i in range(0,n):
                f1 += x[i]**2
            f1 = -0.2 * sqrt(f1/n)
            return f1
        def evalg1(x):
            g1 = zeros(n)
            for i in range(0,n):
                g1[i] = 0.04 * x[i] / (n * evalf1(x))
            return g1

        def evalf2(x):
            f2 = 0
            for i in range(0,n):
                f2 += math.cos(2 * math.pi * x[i])
            f2 = f2/n
            return f2
        def evalg2(x):
            g2 = zeros(n)
            for i in range(0,n):
                g2[i] = -2 * math.pi * math.sin(2 * math.pi * x[i]) / n
            return g2

        def evalf(x):
            f = -20 * math.exp(evalf1(x)) - math.exp(evalf2(x)) + 20 + e
            return f
        def evalg(x):
            g = -20 * math.exp(evalf1(x)) * evalg1(x) - math.exp(evalf2(x)) * evalg2(x)
            return g

    else:
        print('Error!')
        return 0

    xi = 1e-4 #アルミホ条件で用いるξ
    rho = 0.5 #バックトラック法で用いるρ
    epsilon = n * 1e-6 #終了条件で用いるε
    k = 0 #反復回数
    f_num = 0 #目的関数の評価回数(evalfの実行回数)
    memo = zeros(3) #勾配や降下方向のノルムの記録用変数

    if problem == 3:
        iteration = zeros(5) #5回の反復回数を記録する配列
        for i in range(0,5): #行列Aを変えて5回計算を繰り返す
            k = 0

            Z = random.rand(n,n)
            A = dot(Z.T,Z) #正定値行列Aの生成
            x = ones(n)

            while True:
                if k >= 200000: #反復回数が20万回を超えたらbreak
                    print('k =',k,'\n許容最大反復回数を超えました.')
                    return 0

                if method == 1: #最急降下法
                    g = evalg(x)
                    if linalg.norm(g) <= epsilon: #勾配のノルムがε以下になったらbreak
                        break
                    d = - g

                elif method == 2: #ニュートン法
                    g = evalg(x)
                    if linalg.norm(g) <= epsilon: #勾配のノルムがε以下になったらbreak
                        break
                    tau = 0
                    B = evalh(x)
                    I = identity(n)
                    while True:
                        try: #Lをコレスキー分解し降下方向を計算
                            L = linalg.cho_factor(B + tau * I)
                            d = linalg.cho_solve(L,-g)
                            break
                        except: #Lがコレスキー分解出来ないときは単位行列のτ倍を足して正定値にする
                            if tau == 0:
                                tau = 2
                            else:
                                tau = tau * 2

                else:
                    print('Error!')
                    return 0

                t = 1 #ステップサイズの変数
                f = evalf(x)
                dg = inner(d,evalg(x))
                while True: #ステップサイズの決定
                    f_step = evalf(x + t * d)
                    if f_step <= f + xi * t * dg: #アルミホ条件を満たしたらbreak
                        break
                    if search == 1: #ρを掛けてtを更新
                        t = rho * t
                    elif search == 2: #二次補間法を用いてtを更新
                        t_prev = t
                        t = - dg * t**2 / (2 * (f_step - f - dg * t))
                        if t < 0.1 * t_prev or t > 0.9 * t_prev:
                            t = t_prev / 2.0
                    else:
                        print('Error!')
                        return 0
                x = x + t * d #xの更新
                k = k + 1 #反復回数を1増やす
            iteration[i] = k #反復回数を記録
        med = median(iteration) #5回の反復回数の中央値を計算

    else:
        while True:
            if k >= 200000: #反復回数が20万回を超えたらbreak
                print('k =',k,'\n許容最大反復回数を超えました.')
                return 0

            g = evalg(x)
            if method == 1: #最急降下法
                if problem == 4 or problem == 5:
                    print('Error!')
                    return 0
                if linalg.norm(g) <= epsilon: #勾配のノルムがε以下になったらbreak
                    break
                d = - g

            elif method == 2: #ニュートン法
                if problem == 4 or problem == 5:
                    print('Error!')
                    return 0
                if linalg.norm(g) <= epsilon: #勾配のノルムがε以下になったらbreak
                    break
                tau = 0
                B = evalh(x)
                I = identity(n)
                while True:
                    try: #Lをコレスキー分解し降下方向を計算
                        L = linalg.cho_factor(B + tau * I)
                        d = linalg.cho_solve(L,-g)
                        break
                    except: #Lがコレスキー分解出来ないときは単位行列のτ倍を足して正定値にする
                        if tau == 0:
                            tau = 2
                        else:
                            tau = tau * 2

            elif method == 3: #射影勾配法
                if problem == 1 or problem == 2:
                    print('Error!')
                    return 0
                s = 1
                d = fmax(l,fmin(x - s * g,u)) - x

                norm_d = linalg.norm(d)
                memo[0] = memo[1]
                memo[1] = memo[2]
                memo[2] = norm_d #直前3回の降下方向のノルムを記録
                if norm_d <= epsilon: #降下方向のノルムがε以下になったらbreak
                    break

            else:
                print('Error!')
                return 0

            t = 1 #ステップサイズの変数
            f = evalf(x)
            f_num = f_num + 1
            dg = inner(d,g)
            while True: #ステップサイズの決定
                f_step = evalf(x + t * d)
                f_num = f_num + 1
                if f_step <= f + xi * t * dg: #アルミホ条件を満たしたらbreak
                    break
                if search == 1: #ρを掛けてtを更新
                    t = rho * t
                elif search == 2: #二次補間法を用いてtを更新
                    t_prev = t
                    t = - dg * t**2 / (2 * (f_step - f - dg * t))
                    if t < 0.1 * t_prev or t > 0.9 * t_prev:
                        t = t_prev / 2.0
                else:
                    print('Error!')
                    return 0

            if method == 3: #射影勾配法で点列が十分更新されない場合break
                if t * max(abs(d)) <= 1e-16 * max(1,max(abs(x))):
                    break
            x = x + t * d #xの更新
            k = k + 1 #反復回数を1増やす
            if problem == 1 or problem == 2:
                memo[0] = memo[1]
                memo[1] = memo[2]
                memo[2] = linalg.norm(evalg(x)) #直前3回の勾配のノルムを記録

    #用いた手法を出力
    if method == 1:
        print('手法:最急降下法')
    elif method == 2:
        print('手法:ニュートン法')
    else:
        print('手法:射影勾配法')

    #用いた直線探索を出力
    if search == 1:
        print('直線探索:ρ=0.5\n')
    else:
        print('直線探索:2次補間法\n')

    #各課題において実験結果を出力
    if problem == 1:
        print('初期点 x =',x_0,'\n最後の3回の勾配のノルムの値 =',memo)
    if problem == 2:
        print('最後の3回の勾配のノルムの値 =',memo)
    if problem == 3:
        print('次元数 =',n,'\n反復回数の中央値 =',med)
    if problem == 4:
        print('最後の3回の降下方向のノルムの値 =',memo)
    if problem == 5:
        print('初期点 x =',x_0,'\n次元数 =',n,'\n最後の3回の降下方向のノルムの値 =',memo)

    if not problem == 3:
        print('反復回数 =',k,'\n目的関数の評価回数 =',f_num,'\n最適値 f =',evalf(x),'\n解 x =',x)
if __name__ == '__main__': main()
