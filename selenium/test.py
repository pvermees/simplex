import time, os
from selenium import webdriver
from selenium.webdriver.common.keys import Keys

# start test.R first
def boot():
    global d
    d = webdriver.Firefox()
    d.get("http://127.0.0.1:8000")
    
def setup():
    d.find_element('id',"setup").click()

def instrument(inst):
    instrument = d.find_element("id","instrument")
    instrument.clear()
    instrument.send_keys(inst)

def blank(val):
    blank = d.find_element("id","blank")
    blank.clear()
    blank.send_keys(val)
    
def upload(fname):
    upload = d.find_element("id","upload")
    upload.send_keys(fname)

def SHRIMPtest():
    setup()
    instrument("SHRIMP")
    blank("bkg")
    upload("/home/pvermees/Documents/Programming/R/simplex/inst/SHRIMP.op")
