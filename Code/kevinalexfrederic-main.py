import codecs
from sympy import simplify, latex, pprint, pretty, pi
from kevinalexfrederic.models import *
from kevinalexfrederic.plotting import *

# Eenheden:
#   afstanden:  [x] = kpc
#   snelheden:  [v] = km/s
#   massa's:    [m] = M_sol
#   tijden:     [t] = sec . km/kpc (approx. 950e6 jaar)


orbits = [
    # file,    title,                   model,           plot,          plot_opt
    ["box",    "Box",                   SimpleGalaxy,    OrbitPlot,     {},
     #v0, a,   b,    int_period, int_steps, x0, y0, z0, vx0, vy0, vz0, mod_opt
     220, 0.8, 0.75,  10.0,      10000,      2,  2,  2,   0,   0,  0,  {}],

    # file,    title,                   model,           plot,          plot_opt
    ["mass",   "Massaverdeling",        SimpleGalaxy,    DensityPlot,   {},
     #v0, a,   b,    int_period, int_steps, x0, y0, z0, vx0, vy0, vz0, mod_opt
     220, 0.8, 0.75,  10.0,      10000,      2,  2,  2,   0,   0,   0, {}],

    # file,    title,                   model,           plot,          plot_opt
    ["long",   "Long Axis",             SimpleGalaxy,    OrbitPlot,     {},
     #v0, a,   b,    int_period, int_steps, x0, y0, z0, vx0, vy0, vz0, mod_opt
     220, 0.8, 0.69, 200.0,      10000,     20, 50, 40,  10, -75, 100, {}],

    # file,    title,                   model,           plot,          plot_opt
    ["short",  "Short Axis",            SimpleGalaxy,    OrbitPlot,     {},
     #v0, a,   b,    int_period, int_steps, x0, y0, z0, vx0, vy0, vz0, mod_opt
     220, 0.8, 0.6,  200.0,      10000,     20, 50, 40, 100, -75,  10, {}],

    # file,    title,                   model,           plot,          plot_opt
    ["box_bh", "Box met black hole",    BlackHoleGalaxy, BlackHolePlot, {'masses': [5e8/16, 5e8/8, 5e8/4, 5e8/2, 5e8, 10e8]},
     #v0, a,   b,    int_period, int_steps, x0, y0, z0, vx0, vy0, vz0, mod_opt
     220, 0.8, 0.75,  10.0,      10000,      2,  2,  2,   0,   0,   0, {'e': 0.001, 'm': 5e8}],

    # file,    title,                   model,           plot,          plot_opt
    ["bh",     "Box met black hole $M_{\\bullet}$ = 10e9 $M_{\odot}$",    BlackHoleGalaxy, OrbitPlot,     {},
     #v0, a,   b,    int_period, int_steps, x0, y0, z0, vx0, vy0, vz0, mod_opt
     220, 0.8, 0.75, 10.0,      10000,      2,  2,  2,   0,   0,   0, {'e': 0.001, 'm': 10e9}]
]

print "Dichtheid afleiden"

x, y, z, v0, a, b, x0, y0, z0, vx0, vy0, vz0, G = symbols('x y z v0 a b x0 y0 z0 vx0 vy0 vz0 G')
model = SimpleGalaxy(v0, a, b, x0, y0, z0, vx0, vy0, v0)
p = model._potexpr(x, y, z)
rho = 1 / (4 * pi * G) * simplify(diff(p, x, 2) + diff(p, y, 2) + diff(p, z, 2))

print "Uitdrukking voor massadichtheid:"

print
pprint(rho, wrap_line=False)
print

print "Uitdrukking ook opgeslagen in \"dichtheid-uitdrukking.txt\" en \"dichtheid-uitdrukking.tex\""
f = codecs.open("dichtheid-uitdrukking.txt", "w", "utf-8")
f.write(pretty(rho, wrap_line=False, use_unicode=True))
f.close()

f = open("dichtheid-uitdrukking.tex", "w")
f.write(latex(rho))
f.close()

print

print "Plots maken"

num = 1
for orbit in orbits:
    filename, title, model, plotter, plot_opt, v0, a, b, period, steps, x0, y0, z0, vx0, vy0, vz0, mod_opt = orbit

    print "Plot " + str(num) + " van " + str(len(orbits)) + ": \"" + title + "\", " + str(model) + ", " + str(plotter)

    i = model(v0, a, b, x0, y0, z0, vx0, vy0, vz0, **mod_opt)

    plot = plotter(i, **plot_opt)
    plot.plot(title, period, steps)
    print "\tOpslaan als \"" + filename + ".png\""
    plot.save(filename + ".png")
    plot.close()

    num += 1
    print  # lege lijn

print "Klaar :-)"