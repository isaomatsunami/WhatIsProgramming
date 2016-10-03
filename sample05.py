#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
class Triangle:
  def __init__(self,a,b,c):
    self.a = a
    self.b = b
    self.c = c
  def getArea(self):
    # ヘロンの公式
    s = (self.a + self.b + self.c) / 2.0
    return math.sqrt( s * (s - self.a) * (s - self.b) * (s - self.c) )
class Rectangle:
  def __init__(self,a,b):
    self.a = a
    self.b = b
  def getArea(self):
    return self.a * self.b

obj1 = Triangle(3,4,5)
obj2 = Rectangle(6,5)
print( "obj1の面積は" + str( obj1.getArea() ) )
print( "obj2の面積は" + str( obj2.getArea() ) )
