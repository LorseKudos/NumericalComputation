#!/usr/python
# coding: UTF-8
# モンテカルロ法による円周率の計算

from random import *
from math   import *
# 計算結果比較用に組み込み関数により定めた円周率
pi_0 = 2*acos( 0 )

seed()        # 乱数の種の初期化
im = 220  # 乱数発生回数
n = 0.0       # 半径１未満の件数
for i in range(im):
    x = random()
    y = random()
    if ( x**2 + y**2 < 1 ):
      n += 1
pi = 4*n / im

# 結果出力
print 'pi   = ' + str( pi )         # 計算結果
print 'pi_0 = ' + str( pi_0 )       # 組み込み関数の円周率
print 'diff = ' + str( pi - pi_0 )  # 誤差