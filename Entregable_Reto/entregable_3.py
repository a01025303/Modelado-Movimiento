import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from PIL import Image

def no_drag(v0, grados, y0, t, g = 9.81):
    #Convertir ángulo a radianes
    rad = np.radians(grados)
    
    #Velocidad inicial
    v0x = v0 * np.cos(rad)
    v0y = v0 * np.sin(rad)

    #r(t)
    rx = v0x * t
    ry = y0 + v0y * t - (g/2) * t ** 2

    return rx , ry


def coefficient(rho, c, a):
    d = (rho * c * a)/2
    return d

def mass(r, densidad): 
    vol = (4/3) * (np.pi) * r ** 3
    m = vol * densidad
    return m

def maximos(grados, v0):
    rad = np.radians(grados)
    v0y = v0*np.sin(rad)
    v0x = v0*np.cos(rad)
    c = [-4.9, v0y, y0]
    root = np.roots(c)
    for i in range(1):
        maxx= v0x*root[i]
    tmy = v0y/9.8
    maxy = -4.9*tmy**2 + v0y*tmy + y0
    return maxx, maxy

def puntos_c(rx, ry):
    max_list = []
    m0 = (ry[1]-ry[0])/(rx[1]-rx[0])
    for i in range(2, len(ry)):
        m = (ry[i]-ry[i-1])/(rx[i]-rx[i-1])
        if (m0 > 0 and m < 0) or (m0 < 0 and m > 0):
            max_list.append(i)
        m0 = m
    return max_list

def max_drag(x_list, y_list):
    ydmax = max(y_list)
    win = 0
    winner = 0
    for e in range (len(y_list)):
        ind = y_list[e]
        if ydmax == ind:
            win = e
        if ind < 0:
            winner = e-1
            break
    dragxm = x_list[win]
    dragym = y_list[win]
    xdmax = x_list[winner]
    return dragxm, dragym, xdmax

def drag(v0, grados, y0, m, d, dt, itera):
    g = 9.81
    #Grados a radianes
    rad = np.radians(grados)

    # velocidad inicial
    v0x = v0 * np.cos(rad)
    v0y = v0 * np.sin(rad)

    #tiempo
    t = 0
    t_list = [t]

    #velocidad
    v = v0
    vx = v0x
    vy = v0y

    v_list = [v]
    v_x_list = [v0x]
    v_y_list = [v0y]

    # Posición
    x = 0
    y = y0

    x_list = [x]
    y_list = [y]

    # Aceleración inicial
    ax = -(d/m)*v*vx
    ay = -g-(d/m)*v*vy

    a_x_list = [ax]
    a_y_list = [ay]
    
    for _ in range(itera):
            
        # Velocidades para 'x' y 'y' del siguiente paso
        v_x_next = vx + (ax) * dt
        v_y_next = vy + (ay) * dt
        
        # Magnitud del vector de velocidad con componentes v_x_next y v_y_next
        v_next = np.sqrt((v_x_next) ** 2 + (v_y_next) ** 2)

        # Agregar valores a listas v_list, v_x_list, v_y_list. 
        v_list.append(v_next)
        v_x_list.append(v_x_next)
        v_y_list.append(v_y_next)

        # Posiciones 'x' y 'y' del siguiente paso
        x_next = x + v_x_next * dt + (1/2) * ax * (dt ** 2)
        y_next = y + v_y_next * dt + (1/2) * ay * (dt ** 2)
        
        # Agregar calores a listas x_list y y_list
        x_list.append(x_next)
        y_list.append(y_next)
        
        # Aceleraciones para 'x' y 'y' del siguiente paso
        a_x_next = -(d/m) * v * v_x_next
        a_y_next = -g -(d/m) * v * v_y_next

        # Agregar valores a listas a_x_list y a_y_list
        a_x_list.append(a_x_next)
        a_y_list.append(a_y_next)
        
        vx = v_x_next
        vy = v_y_next
        v = v_next

        x = x_next
        y = y_next

        ax = a_x_next
        ay = a_y_next

        # Calcular tiempo y guardarlo en una lista t_list
        t += dt
        t_list.append(t)
        
    return x_list, y_list, v_list, v_x_list, v_y_list, a_x_list, a_y_list, t_list

def trayectorias1(rx, ry, x_list, y_list, maxx, maxy, col):
    fig = plt.figure(figsize = (15, 9))
    plt.plot(rx, ry, label = "Sin resistencia al aire")
    c = ""
    if col == "Morado":
        c = "purple"
    if col == "Verde":
        c = "green"
    if col == "Negro":
        c = "black"
    if col == "Rojo":
        c = "red"
    plt.plot(x_list, y_list, label = "Con resistencia al aire", color = c)
    plt.grid()
    plt.axis("scaled")
    plt.title("Trayectorias con y sin resistencia al aire")
    plt.xlim([0, maxx + 500])
    plt.ylim([0, maxy + 500])
    _ = plt.legend()
    return fig

def tray_extra(t_list, v_x_list, v_y_list, v_list, a_x_list, a_y_list, check):
    figura = plt.figure(figsize= (15, 8))
    if check == "Velocidad vs. tiempo":
        plt.title("Velocidad contra tiempo")
        plt.plot(t_list, v_x_list, label="Componente x")
        plt.plot(t_list, v_y_list, label="Componente y", color = "green")
        plt.plot(t_list, v_list, label="Velocidad", color = "purple")
    if check == "Aceleración vs. tiempo":
        plt.title("Aceleración contra tiempo")
        plt.plot(t_list, a_x_list, label="Componente x")
        plt.plot(t_list, a_y_list, label="Componente y", color = "green")
    plt.xlim([0, 100])
    plt.grid()
    _ = plt.legend()
    return figura


#Simulación Streamlit
st.title("Simulación de erupción volcánica")
st.subheader("Ana Paula Katsuda, A01025303")
st.write("¡Binvenid@ a la simulación de erupción volcánica del Popocatépetl! Cambia los valores para visualizar la trayectoria de los proyectiles con y sin resistencia al aire.")
st.write("Nota: la trayectoria sin resistencia al aire está marcada con azul.")
foto = Image.open('popo.jpg')
st.image(foto, caption="Popocatépetl (Chacón, L., 2019)")
with st.beta_expander("Más Información"):
    st.write("El Popocatépetl está ubicado a 70 km de la Ciudad de México a una altura de 5452 m sobre el nivel del mar. El volcán es considerado uno de los más peligrosos de México dado su historial de erupciones. Asimismo, implica un riesgo puesto a que un aproximado de 500,000 personas viven a 10-30 km de distancia del mismo (Ibargüengoitia, M., 2011).")
    st.write("Actualmente, las erupciones más recientes del Popocatépetl son consideradas como moderadas. El volcán estuvo inactivo durante 70 años y reinició su actividad en 1994. Desde entonces, sus cenizas han logrado alcanzar Puebla,la Ciudad de México y otras poblaciones (Osorio, M., Puente, L., Valdés, C., 2014)")
with st.beta_expander("¿De dónde salieron los datos?"):
    st.write("Se ha reportado que los balísticos del Popocatépetl son de 0.20 a 0.60 m de diámetro, y de 2100-2600 kg/m^3. Además, la velocidad inicial de los mismos se encuentra entre 180 y 230 m/s (Taddeucci, J., et al, 2017) y el ángulo incial típico es de 25 a 70º (Rodríguez, D., Córdoba, G. y Costa, A., 2018). Por otro lado, se han registado densidades de aire desde 0.1 (Yunus, Ç., Cimbala, J., 2006) hasta 1.2 kg/m^3 (Pacheco, J., 2008) y un coeficiente de arrastre desde 0.1 hasta 1.25 (Taddeucci, J., et al, 2017).")

st.write("A continuación se muestra la trayectoria de los proyectiles tanto con resistencia al aire como sin resistencia al aire. Asimismo, es posible visualizar la altura máxima a la que llega cada proyectil y sus coordenadas.")

#Datos proyecil sin resistencia al aire
st.sidebar.subheader("Datos iniciales para proyectil sin resistencia al aire")
v0 = st.sidebar.slider ("Velocidad inicial (m/s)", 180, 230)
grados = st.sidebar.slider("Grados º", 25, 80)
y0 = st.sidebar.slider("Altura inicial (m)", 3000, 5500, 3315)
t = np.linspace(0, 70, num = 60)
#Datos proyectil con resistencia al aire
st.sidebar.subheader("Datos iniciales para proyectil con resistencia al aire")
r = st.sidebar.number_input("Radio del proyectil (m)", min_value=0.1, max_value=0.5)
a = np.pi*r**2
rho = st.sidebar.number_input("Densidad del aire", min_value= 0.1, max_value=1.3)
densidad = st.sidebar.slider("Densidad del proyectil", min_value=2100, max_value=2600)
c = st.sidebar.number_input("Coeficiente de arrastre", min_value = 0.1, max_value=1.25)
itera = 10000
dt = 0.1
colores = ["Morado", "Negro", "Verde", "Rojo"]
col = st.sidebar.selectbox("Color de trayectoria con resistencia al aire", colores)

#Llamada de funciones
rx, ry = no_drag(v0, grados, y0, t, g=9.81)
m = mass(r, densidad)
d = coefficient(rho, c, a)
[x_list, y_list, v_list, v_x_list, 
 v_y_list, a_x_list, a_y_list, t_list] = drag(v0, grados, y0, m, d, dt, itera)
maxx, maxy = maximos(grados, v0)
max_list = puntos_c(rx, ry)
dragxm, dragym, xdmax = max_drag(x_list, y_list)
fig = trayectorias1(rx, ry, x_list, y_list, maxx, maxy, col)

#Señalar punto máximo
for indx in max_list:
    plt.scatter(rx[indx-1],ry[indx-1],c='b')
plt.scatter(dragxm, dragym,c='black')
#Imprimir primer gráfico
st.pyplot(fig)
with st.beta_expander("Ver datos"):
    st.write("Alcance máximo de proyectil sin resistencia al aire: ", maxx)
    st.write("Altura máxima de proyectil sin resistencia al aire: ", maxy)
    st.write("Alcance máximo de proyectil con resistencia al aire: ", xdmax)
    st.write("Altura máxima de proyectil con resistencia al aire: ", dragym)
st.write("Es posible notar una diferencia considerable en ambas trayectorias. La trayectoria con resistencia al aire tiene una altura máxima y un alcance máximo menor a la trayectoria sin resistencia. Lo anterior debido a que la resistencia al aire se opone al movimiento.")
st.header("Gráficos adicionales para proyectil con resistencia al aire")
#Opciones siguiente gráfico
opciones = ["Velocidad vs. tiempo", "Aceleración vs. tiempo"]
check = st.radio("Elige el gráfico: ", opciones)
figura = tray_extra(t_list, v_x_list, v_y_list, v_list, a_x_list, a_y_list, check)

#Imprimir gráfico
st.pyplot(figura)

#Referencias
with st.beta_expander("Referencias"):
    st.write("Alatorre, M. (2011). A model of volcanic explosions at Popocatépetl volcano (Mexico)... [Sitio Web]. Recuperado de: https://d-nb.info/1015203132/34")
    st.write("Becerra, L., Guardado, M. (2001). Estimación de la incertidumbre en la determinación de la densidad del aire. CENAM. [Sitio Web]. Recuperado de: http://www.cenam.mx/myd/DENSIDAD%20DEL%20AIRE%20abril-20031.pdf")
    st.write("Bucheli, E. Modelación Computacional del Movimiento. Tecnológico de Monterrey, Ciudad de México. 7 oct. 2020")
    st.write("Chacón, L. (2019). Popocatépetl: gigante bajo vigilancia. [Imagen]. Recuperado de:  https://www.infochannel.info/popocatepetl-gigante-bajo-vigilancia")
    st.write("INEGI (2017). Anuario estadístico y geográfico de Puebla 2017. [Sitio Web]. Recuperado de: https://www.datatur.sectur.gob.mx/ITxEF_Docs/PUE_ANUARIO_PDF.pdf")
    st.write("Osorio, M., Puente, L., Valdés, C. (2014). Historia de la actividad del Volcán Popocatépetl. [Sitio Web]. Recuperado de: http://www.cenapred.gob.mx/es/Publicaciones/archivos/225-HISTORIADELAACTIVIDADDELVOLCNPOPOCATPETL-17AOSDEERUPCIONES.PDF")
    st.write("Pacheco, J. (2008) Análisis y Comparación de las emisiones… UNAM. [Sitio Web]. Recuperado de: http://www.ptolomeo.unam.mx:8080/xmlui/bitstream/handle/132.248.52.100/8070/Tesis_Completa.pdf?sequence=1")
    st.write("Shiffman, D. (s.f.). Resistencia del aire y de fluidos. Khan Academy. [Sitio Web]. Recuperado de: https://es.khanacademy.org/computing/computer-programming/programming-natural-simulations/programming-forces/a/air-and-fluid-resistance")
    st.write("Taddeucci, J., Alatorre‐Ibargüengoitia, M. A., Cruz‐Vázquez, O., Del Bello, E., Scarlato, P., and Ricci, T. (2017), In-flight dynamics of volcanic ballistic projectiles. [Sitio Web]. Recuperado de: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017RG000564")
    st.write("Yunus, Ç., Cimbala, J. (2006). Mecánica de fluidos: Fundamentos y aplicaciones. [Sitio Web]. Recuperado de: http://tesis.pucp.edu.pe/repositorio/bitstream/handle/20.500.12404/5421/SOTOMAYOR_DENIS_SIMULACION_NUMERICA_INTERCAMBIADOR_CALOR_FLUJO_TRANSVERSAL_ALETEADO_ANEXOS.pdf?sequence=2&isAllowed=y")