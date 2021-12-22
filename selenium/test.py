def start():
    import time, os
    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys
    driver = webdriver.Firefox()
    os.system("Rscript test.R &")
    driver.get("http://127.0.0.1:8000")

def refresh():
    driver.refresh()

def setup():
    driver.find_element('id',"setup").click()

def instrument(inst):
    instrument = driver.find_element("id","instrument")
    instrument.clear()
    instrument.send_keys(inst)

def blank(val):
    blank = driver.find_element("id","blank")
    blank.clear()
    blank.send_keys(val)
    
def upload(fname):
    upload = driver.find_element("id","upload")
    upload.send_keys(fname)

def SHRIMPtest():
    setup()
    instrument("SHRIMP")
    blank("bkg")
    upload("/home/pvermees/Documents/Programming/R/simplex/inst/SHRIMP.op")

SHRIMPtest()
