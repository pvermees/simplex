import time, os, glob
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select

# start test.R first, then start Python (Ctrl-c, Ctrl-p) and load script (Ctrl-l)
# to start test.py: boot(), to refresh: d.refresh()
def boot():
    global d
    d = webdriver.Firefox()
    d.get("http://127.0.0.1:8000")

def click(identifier):
    d.find_element('id',identifier).click()

def select(parent,daughter):
    Select(d.find_element("id",parent)).select_by_value(daughter)

def enter(identifier,val):
    item = d.find_element("id",identifier)
    item.clear()
    item.send_keys(val)
    item.send_keys(Keys.ENTER)
    
def upload(fname):
    upload = d.find_element("id","upload")
    upload.send_keys(fname)

def mupload(wildcard):
    fnames = glob.glob(wildcard)
    unames = fnames[0]
    for i in range(1,len(fnames)):
        unames += "\n" + fnames[i]
    d.implicitly_wait(10)
    upload(unames)

def SHRIMPtest():
    click("setup")
    enter("instrument","SHRIMP")
    enter("ions","Zr2O,Pb204,bkg,Pb206,Pb207,Pb208,U238,ThO,UO,UO2")
    enter("bkg","bkg")
    upload("/home/pvermees/Documents/Programming/R/simplex/inst/SHRIMP.op")
    time.sleep(1)
    click("drift")

def UPbtest():
    click("setup")
    mupload("/home/pvermees/Documents/SIMS/20210528 Tanz U-Th-Pb/*.asc")
    enter("ions","90Zr2O,92Zr2O,200.5,94Zr2O,Pb204,Pb206,Pb207,Pb208,HfO2,Th232,U238,ThO,UO,UO2")
    time.sleep(1)
    click("drift")
    time.sleep(5)
    click("logratios")

def Otest():
    click("setup")
    mupload("/home/pvermees/Documents/SIMS/20200402 Tanghejun L194 O test/*.asc")
    enter("ions","O16,17,O18")
    time.sleep(1)
    enter("num","O18")
    time.sleep(1)
    enter("den","O16")
    time.sleep(1)
    click("drift")
    time.sleep(1)
    click("logratios")
    time.sleep(1)
    click("calibration")
    
def geochroncalibrationtest():
    click("calibration")
    select("standtype","measured")
    select("standcomp","t2stand")
    t = d.find_element("id","t")
    t.clear()
    t.send_keys(500)
    st = d.find_element("id","st")
    st.clear()
    st.send_keys(10)
    click("convert-t2stand")

def oxygencalibrationtest():
    click("setup")
    select("methods","IGG-O")
    enter("num","O18")
    enter("den","O16")
    click("drift")
    time.sleep(1)
    click("logratios")
    time.sleep(1)
    click("calibration")
    time.sleep(1)
    click("calibrate")
    time.sleep(1)
    click("samples")

def oxygenfinishtest():
    click("setup")
    select("methods","IGG-O")
    time.sleep(1)
    click("calibration")
    time.sleep(1)
    click("calibrate")
    time.sleep(1)
    click("samples")
    time.sleep(1)
    click("calibrate")

def PbPbtest():
    click("setup")
    enter("num","Pb204,Pb207,Pb208")
    time.sleep(1)
    enter("den","Pb206,Pb206,Pb206")
    time.sleep(1)
    click("drift")
    time.sleep(3)
    click("logratios")
    time.sleep(2)
    click("calibration")
