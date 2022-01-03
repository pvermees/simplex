import time, os, glob
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select

# start test.R first
def boot():
    global d
    d = webdriver.Firefox()
    d.get("http://127.0.0.1:8000")

def click(id):
    d.find_element('id',id).click()

def select(parent,daughter):
    Select(d.find_element("id",parent)).select_by_value(daughter)
    
def instrument(inst):
    instrument = d.find_element("id","instrument")
    instrument.clear()
    instrument.send_keys(inst)

def blank(val):
    blank = d.find_element("id","blank")
    blank.clear()
    blank.send_keys(val)
    blank.send_keys(Keys.ENTER);

def ions(val,identifier="ions"):
    ions = d.find_element("id",identifier)
    ions.clear()
    ions.send_keys(val)
    ions.send_keys(Keys.ENTER);
    
def upload(fname):
    upload = d.find_element("id","upload")
    upload.send_keys(fname)

def SHRIMPtest():
    click("setup")
    instrument("SHRIMP")
    blank("bkg")
    upload("/home/pvermees/Documents/Programming/R/simplex/inst/SHRIMP.op")
    click("drift")

def Camecatest():
    click("setup")
    fnames = glob.glob("/home/pvermees/Documents/SIMS/20210528 Tanz U-Th-Pb/*.asc")
    unames = fnames[0]
    for i in range(1,len(fnames)):
        unames += "\n" + fnames[i]
    d.implicitly_wait(10)
    upload(unames)
    ions("90Zr2O,92Zr2O,200.5,94Zr2O,Pb204,Pb206,Pb207,Pb208,HfO2,Th232,U238,ThO,UO,UO2")
    time.sleep(1)
    click("drift")
    time.sleep(5)
    click("logratios")
    
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
    click("convert")

def oxygencalibrationtest():
    click("setup")
    select("methods","IGG-O")
    ions("O18",identifier="num")
    ions("O16",identifier="den")
    click("drift")
    time.sleep(0.2)
    click("logratios")
    time.sleep(0.2)
    click("calibration")
    select("standcomp","del2stand")
    
