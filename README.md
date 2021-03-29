# Risk-and-Misallocation

Note this is in Markdown:  https://www.markdownguide.org/basic-syntax/

## High level solution method notes & questions 

For solving I think we are basically looking at a double loop, "guess w, guess r, solve PE problem, find the ergodic distribution, update r, update w". This double loop makes it significantly more annoying than a regular HA model which usually only has one loop. Even the model I worked with before essentially only had one loop... because there were positive gov't bonds and this wasn't important to dynamics, it was easiest to just pick r and assume the gov't bond supply adjusted to the rate). I think the model you solved before was also a double loop, so we could just use that? 

But in any case I think probably worthwhile to spend a bit of time up front thinking about the methods we use, since I anticipate resolving this or more complicated variations a lot. I think there are a few major choices to make: 

1. I see you used bisection to update r & w. Not sure if there is a better way to do this. 
2. I see you used VFI to solve the PE problem. I think there are lots of ways to do this and VFI is one of the slower ways? Have you seen these notebooks Chase wrote on QuantEcon? I would read the Julia or the Python one since Chase complained that Matlab didn't have Jupyter integration and also I know Chase hated Matlab, but I imagine the relative algorithm speeds shouldn't vary too much across programs. (There are also papers that do this sort of exercise but I find Chase's notebooks very practical... did you meet Chase? Nice guy. Very prosocial. Didn't really care about economics but cared about coding...) 
    - For Moll, I don't think there is a double loop of this type (a risk wedge) up on his site, although there is HANK sort of stuff (ie an illiquid asset). 
    - Want to look around a little more 
3. Updating the ergodic distribution is also often a choke point, although I am not so sure about different methods to do this. 

