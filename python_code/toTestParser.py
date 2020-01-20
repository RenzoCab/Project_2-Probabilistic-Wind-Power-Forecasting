#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse


# In[2]:


parser = argparse.ArgumentParser(description='Likelihood Evaluator v1.0')
parser.add_argument('-filename', help=' Config file name or path')
parser.add_argument('--version', action='version', version='Likelihood Evaluator v1.0')
args = parser.parse_args()

print(args)


# In[ ]:




