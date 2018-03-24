#!/use/bin/env python
# -*- coding: UTF-8 -*-

from numpy import *
from scipy import *
from scipy import linalg
import string
import sys

def main():
    args = sys.argv[1:]
    search = string.atoi(args[0]) #直線探索を決める変数

    m = int(raw_input('反復回数を出力したいときは1 評価回数を出力したいときは2を入力してください :'))
    for i in range(0,100):
        n = i+1
        x = ones(n) / 3

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

        xi = 1e-4 #アルミホ条件で用いるξ
        rho = 0.5 #バックトラック法で用いるρ
        epsilon = n * 1e-6 #終了条件で用いるε
        k = 0 #反復回数
        f_num = 0 #目的関数の評価回数(evalfの実行回数)

        while True:
            if k >= 200000: #反復回数が20万回を超えたらbreak
                print 'k =',k,'\n許容最大反復回数を超えました.'
                return 0

            g = evalg(x)

            s = 1
            d = fmax(l,fmin(x - s * g,u)) - x

            norm_d = linalg.norm(d)
            if norm_d <= epsilon: #降下方向のノルムがε以下になったらbreak
                break

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

            if t * max(abs(d)) <= 1e-16 * max(1,max(abs(x))):
                break
            x = x + t * d #xの更新
            k = k + 1 #反復回数を1増やす
        if m == 1:
            print n,k
        elif m == 2:
            print n,f_num

if __name__ == '__main__': main()
