from ctypes import c_double, c_int

from ROOT import TFile, TTree
import ROOT

guan = TFile("guan.root", "recreate")

ming = TTree("ming", "This is my first tree")

i = c_int(1)
a = c_double(2.2)
b = c_double(0.0)

ming.Branch("i", i, "i/I")  # "i" 为分支的名称，i为地址（c_int本身就是地址）,"i/I"为i的类型
ming.Branch("a", a, "a/D")
for j in range(100):
    i = j
    a = i * 0.1
    b = a * a * 9
    ming.Fill()
ming.Write()

guan.ls()
