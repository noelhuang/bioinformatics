import sys


# Determines the percentage of characters that is 'a' or 'b' or 'c' from a string input
def percentage_a(input_text):
    abc_count = 0
    abc_list = ['a', 'b', 'c']
    for i in input_text:
        if i in abc_list:
            abc_count = abc_count + 1
    abc_content = abc_count / len(input_text)
    return abc_content
