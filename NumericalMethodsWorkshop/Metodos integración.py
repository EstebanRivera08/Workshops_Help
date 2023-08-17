
from numpy import *



def punto_medio_comp(f,a,b,N) :
    h = (b-a)/N
    x = linspace(a,b,N+1)
    x_prom = (x[2:N+1]+x[1:N])/2
    g = h*sum(f(x_prom))
    return g

def simpson_comp(f,a,b,N) :
    h = (b-a)/N
    x = linspace(a,b,N+1)
    Impar = Par = 0
    for i in range(N) :
        if ((-1)^i)==1 :
           Impar=Impar+f(x[i])
        else:
           Par=Par+f(x[i])

    g = h/3*(f(a)+4*Impar+2*Par+f(b))
    return g

def trapecio_comp(f,a,b,N) :
    h = (b-a)/N
    x = linspace(a,b,N+1)
    g = h/2*(f(a)+f(b)+2*sum(f(x[2:N])))
    return g

def integraCuadGauss2p(funcionx,a,b):
    x0 = -1/np.sqrt(3)
    x1 = -x0
    xa = (b+a)/2 + (b-a)/2*(x0)
    xb = (b+a)/2 + (b-a)/2*(x1)
    area = ((b-a)/2)*(funcionx(xa) + funcionx(xb))
    return(area)

# INGRESO
fx = lambda x: np.sqrt(x)*np.sin(x)

# intervalo de integración
a = 1
b = 3
tramos = 4

def Cuadratura_gaussiana(fx,a,b,n) :
# PROCEDIMIENTO
    m = n+1
    xi = np.linspace(a,b,m)
    area = 0
    for i in range(0,n,1):
        deltaA = integraCuadGauss2p(fx,xi[i],xi[i+1])
        area = area + deltaA
# SALIDA
    print('Integral: ', area)

# GRAFICAR por cada Segmento
# para concepto con pocos segmentos
    x0 = - 1 /np.sqrt(3)
    x1 = 1/ np.sqrt(3)

# arregos para gráficas
    xi = np.array([])
    fi = np.array([])

    xat = np.array([])
    xbt = np.array([])

    recta = np.array([])

    mtramo = 10
    subtramo = np.linspace(a, b, m)

    for i in range(0, n, 1):
        at = subtramo[i]
        bt = subtramo[i + 1]

        xit = np.linspace(at, bt, muestrastramo)
        fit = fx(xit)

        xi = np.concatenate((xi, xit))
        fi = np.concatenate((fi, fit))

    # puntos xa y xb por tramo
        xa = (bt + at) / 2 + (bt - at) / 2 * (x0)
        xb = (bt + at) / 2 + (bt - at) / 2 * (x1)

        xat = np.concatenate((xat, [xa]))
        xbt = np.concatenate((xbt, [xb]))

    # Recta entre puntos x0 y x1 por tramo
        pendiente = (fx(xb) - fx(xa)) / (xb - xa)
        b0 = fx(xa) - pendiente * xa
        linea = b0 + pendiente * xit
        recta = np.concatenate((recta, linea))

# Marcadores 'o' de xa y xb por tramos
    puntox = np.concatenate((xat, xbt))
    puntoy = fx(puntox)

# Graficando
# Trazado de lineas
    plt.plot(xi, recta, label='grado 1', color='tab:orange')
    plt.fill_between(xi, 0, recta, color='tab:olive')
    plt.plot(xi, fi, label='f(x)', color='blue')

# Verticales para dividir los tramos
    for j in range(0, len(subtramo), 1):
        plt.axvline(subtramo[j], color='tab:gray')

# Marcadores de puntos xa y xb por tramos
    for j in range(0, len(xat), 1):
        plt.axvline(xat[j], color='w')
        plt.axvline(xbt[j], color='w')

    plt.plot(puntox, puntoy, 'o', color='g')

    plt.title('Integral: Cuadratura Gauss')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()

    plt.show()
