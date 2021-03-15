import mechanize
import re
import os, sys 
import webbrowser

# n.b. HCF = highest common factor and this is the same as the greatest common
# divisor.

def gcd(a, b):
    """Calculate the Greatest Common Divisor of a and b.

    Unless b==0, the result will have the same sign as b (so that when
    b is divided by it, the result comes out positive).
    """
    while b:
        a, b = b, a%b
    return a

# Initialise the browser and set the url of the quiz page.
br = mechanize.Browser()
page = "http://www.quantbet.com/quiz/dev"

# Open the webpage and read the content.
handle = br.open(page)
content = handle.read()
print(content)

# Locate the numbers for the calculation. First all, the numbers on the page
# are extracted and then these are filtered for numbers of more than four
# digits - the two numbers that we need to find the HCF of.

numbers_all = re.findall('\d+', content)
numbers = []  # The numbers required for the HCF calculation
for num in numbers_all:
    if len(num) > 4:
        numbers.append(num)

# Compute the HCF of the two numbers.
hcf = gcd(int(numbers[0]), int(numbers[1]))

# By inspecting the html code of the webpage, it is found that the HCF must be
# entered in a form which has an id of 'quiz' and the input name is 'divisor'.
for form in br.forms():
  if form.attrs['id'] == 'quiz':
    br.form = form
    break

br["divisor"] = str(hcf)
response = br.submit()  # capture the response after submission

# Convert the response into an html string that can be opened in a browser
html = ''
for i in response:
  html += str(i)
print(html)
br.close()

# Open html code in browser
path = os.path.abspath('temp.html')
url = 'file://' + path
with open(path, 'w') as f:
  f.write(html)
webbrowser.open(url)
