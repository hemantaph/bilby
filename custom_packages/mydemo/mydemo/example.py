"""Example module"""
from mydemo.addition import test

def say_hello(name):
    """Says hello to given name"""
    print(f"hello, {name}")
    
def add(a,b):
    return( test.add_(a,b)+square(a) )

def square(i):
    return(i**2)