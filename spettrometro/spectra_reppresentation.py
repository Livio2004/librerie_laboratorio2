import numpy as np

def wave2rgb(wave):
    # This is a port of javascript code from  http://stackoverflow.com/a/14917481
    gamma = 0.8
    intensity_max = 1
 
    if wave < 380:
        red, green, blue = 0, 0, 0
    elif wave < 440:
        red = -(wave - 440) / (440 - 380)
        green, blue = 0, 1
    elif wave < 490:
        red = 0
        green = (wave - 440) / (490 - 440)
        blue = 1
    elif wave < 510:
        red, green = 0, 1
        blue = -(wave - 510) / (510 - 490)
    elif wave < 580:
        red = (wave - 510) / (580 - 510)
        green, blue = 1, 0
    elif wave < 645:
        red = 1
        green = -(wave - 645) / (645 - 580)
        blue = 0
    elif wave <= 780:
        red, green, blue = 1, 0, 0
    else:
        red, green, blue = 0, 0, 0
 
    # let the intensity fall of near the vision limits
    if wave < 380:
        factor = 0
    elif wave < 420:
        factor = 0.3 + 0.7 * (wave - 380) / (420 - 380)
    elif wave < 700:
        factor = 1
    elif wave <= 780:
        factor = 0.3 + 0.7 * (780 - wave) / (780 - 700)
    else:
        factor = 0
 
    def f(c):
        if c == 0:
            return 0
        else:
            return intensity_max * pow (c * factor, gamma)
 
    return f(red), f(green), f(blue)

def draw_spectrum(λ,sigma=None,ytitle=None):
    import matplotlib.pyplot as plt

    height = 1
    width = 10 * height

    plt.style.use('dark_background')  
    fig = plt.figure(figsize=(width, height))
    
    if sigma is not None:
        for l,s in zip(λ,sigma):
            c = wave2rgb(float(l))   
            plt.axvspan(l-s, l+s, color=c, alpha=0.2)
            plt.axvline(x=l, color=c)

    else:
        for l in λ:
            c = wave2rgb(float(l))
            plt.axvline(x=l, color=c)
    
    plt.xticks(np.arange(350, 701, 50))
    plt.yticks([])
    plt.xlabel('λ [nm]')   
    plt.ylabel(ytitle)
  
    plt.show()
    plt.style.use('default')
    
