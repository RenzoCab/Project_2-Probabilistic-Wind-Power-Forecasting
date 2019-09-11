import json
import os
os.chdir('/home/alhaddwt/Insync/waleedhad@gmail.com/Google Drive/GitLab/wind-power/python_code/config')

class Test(object):
    def __init__(self, data):
	    self.__dict__ = json.load(data)
