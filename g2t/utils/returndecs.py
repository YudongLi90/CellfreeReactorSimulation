"""
Author: Yudong Li (Yudong.Li@nrel.gov)
"""
def non_negative_return(func):
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        return max(res, 0)
    return wrapper