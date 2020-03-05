var abs = Math.abs;

function array_fill(i, n, v) {
    var a = [];
    for (; i < n; i++) {
        a.push(v);
    }
    return a;
}

/**
 * Gaussian elimination
 * @param  array A matrix
 * @param  array x vector
 * @return array x solution vector
 */

function gauss(A, x) {

    var i, k, j;

    // Just make a single matrix
    for (i = 0; i < A.length; i++) {
        A[i].push(x[i]);
    }
    var n = A.length;

    for (i = 0; i < n; i++) {
        // Search for maximum in this column
        var maxEl = abs(A[i][i]),
            maxRow = i;
        for (k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }


        // Swap maximum row with current row (column by column)
        for (k = i; k < n + 1; k++) {
            var tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (k = i + 1; k < n; k++) {
            var c = -A[k][i] / A[i][i];
            for (j = i; j < n + 1; j++) {
                if (i === j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    x = array_fill(0, n, 0);
    for (i = n - 1; i > -1; i--) {
        x[i] = A[i][n] / A[i][i];
        for (k = i - 1; k > -1; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }

    return x;
}

class MetodoPolinomial {

    constructor(dX, dY, iteraciones) {
        this.__datosX = dX;
        this.__datosY = dY;
        this.__size = iteraciones;
        // this.__inter = interpolar;
        this.__XX = new Array();
        this.__XXX = new Array();
        this.__XXXX = new Array();
        this.__XY = new Array();
        this.__XXY = new Array();
    }

    get getDataX() {

        return this.__datosX;
    }

    get getDataY() {
        return this.__datosY;
    }

    get getSize() {
        return this.__size;
    }

    get getInter() {
        return this.__inter;
    }

    get getXX() {
        var vectorget = new Array();
        for (let i = 0; i < this.__size; i++) {
            vectorget[i] = this.__datosX[i] * this.__datosX[i];
        }
        this.__XX = vectorget;
        return this.__XX;
    }

    get getXXX() {
        var vectorget = new Array();
        for (let i = 0; i < this.__size; i++) {
            vectorget[i] = this.__datosX[i] * this.__datosX[i] * this.__datosX[i];
        }
        this.__XXX = vectorget;
        return this.__XXX;
    }

    get getXXXX() {
        var vectorget = new Array();
        for (let i = 0; i < this.__size; i++) {
            vectorget[i] = this.__datosX[i] * this.__datosX[i] * this.__datosX[i] * this.__datosX[i];
        }
        this.__XXXX = vectorget;
        return this.__XXXX;
    }

    get getXY() {
        var vectorget = new Array();
        for (let i = 0; i < this.__size; i++) {
            vectorget[i] = this.__datosX[i] * this.__datosY[i];
        }
        this.__XY = vectorget;
        return this.__XY;
    }

    get getXXY() {
        var vectorget = new Array();
        for (let i = 0; i < this.__size; i++) {
            vectorget[i] = this.__datosX[i] * this.__datosX[i] * this.__datosY[i];
        }
        this.__XXY = vectorget;
        return this.__XXY;
    }

    get getSumaX() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i];
        }
        return suma;
    }

    get getSumaY() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosY[i];
        }
        return suma;
    }

    get getSumaXX() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i] * this.__datosX[i];
        }
        return suma;
    }

    get getSumaXXX() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i] * this.__datosX[i] * this.__datosX[i];
        }
        return suma;
    }

    get getSumaXXXX() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i] * this.__datosX[i] * this.__datosX[i] * this.__datosX[i];
        }
        return suma;
    }

    get getSumaXXXXX() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i] * this.__datosX[i] * this.__datosX[i] * this.__datosX[i] * this.__datosX[i];
        }
        return suma;
    }

    get getSumaXY() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i] * this.__datosY[i];
        }
        return suma;
    }

    get getSumaYY() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosY[i] * this.__datosY[i];
        }
        return suma;
    }

    get getSumaXXY() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i] * this.__datosX[i] * this.__datosY[i];
        }
        return suma;
    }

    get getSumaXXXY() {
        var suma = 0;
        for (let i = 0; i < this.__size; i++) {
            suma += this.__datosX[i] * this.__datosX[i] * this.__datosX[i] * this.__datosY[i];
        }
        return suma;
    }

    get imprimirDatos() {
        var datosImprimir = [this.__datosX, this.__datosY, this.getXX, this.getXXX, this.getXXXX, this.getXY, this.getXXY];
        return datosImprimir;
    }

    get getPolinomio() {
        const polinomio = new Array();
        polinomio[0] = [`${this.__size}A0`, `${this.getSumaX}A1`, `${this.getSumaXX}A2`, `=${this.getSumaY}`];
        polinomio[1] = [`${this.getSumaX}A0`, `${this.getSumaXX}A1`, `${this.getSumaXXX}A2`, `=${this.getSumaXY}`];
        polinomio[2] = [`${this.getSumaXX}A0`, `${this.getSumaXXX}A1`, `${this.getSumaXXXX}A2`, `=${this.getSumaXXY}`];
        return polinomio;
    }

    get getMatriz() {
        const polinomio = new Array();
        polinomio[0] = [this.__size, this.getSumaX, this.getSumaXX, this.getSumaY];
        polinomio[1] = [this.getSumaX, this.getSumaXX, this.getSumaXXX, this.getSumaXY];
        polinomio[2] = [this.getSumaXX, this.getSumaXXX, this.getSumaXXXX, this.getSumaXXY];
        return polinomio;
    }

    get getPolinomioResultanteGrado2() {
        const variables = new Array();
        variables[0] = [this.__size, this.getSumaX, this.getSumaXX];
        variables[1] = [this.getSumaX, this.getSumaXX, this.getSumaXXX];
        variables[2] = [this.getSumaXX, this.getSumaXXX, this.getSumaXXXX];
        const Objetivos = [this.getMatriz[0][3], this.getMatriz[1][3], this.getMatriz[2][3]];
        const resultado = gauss(variables, Objetivos);
        return resultado;
    }

    get getSumasRows() {
        let contenedor = new Array();
        contenedor = [this.getSumaX, this.getSumaY, this.getSumaXX, this.getSumaXXX, this.getSumaXXXX, this.getSumaXY, this.getSumaXXY];
        return contenedor;
    }

    get getCorrelacion() {
        let correlacion = 0
        let numerador = 0;
        let denominador = 0;
        numerador = (this.__size * this.getSumaXY) - (this.getSumaX * this.getSumaY);
        var p1 = Math.sqrt(this.__size * this.getSumaXX - this.getSumaX * this.getSumaX);
        var p2 = Math.sqrt(this.__size * this.getSumaYY - this.getSumaY * this.getSumaY);
        denominador = p1 * p2;
        correlacion = numerador / denominador
        return correlacion.toFixed(4);
    }

    get getVariablesABlineal() {
        var a = 0;
        var n1 = this.__size * this.getSumaXY - this.getSumaX * this.getSumaY;
        var b = 0;
        var n2 = this.getSumaY * this.getSumaXX - this.getSumaX * this.getSumaXY;
        var d = this.__size * this.getSumaXX - this.getSumaX * this.getSumaX;
        a = n1 / d;
        b = n2 / d;
        const arrayDatos = [a, b];
        return arrayDatos;
    }

    get getValExp() {
        var c = Math.exp(this.getVariablesABlineal[1]);
        var a = this.getVariablesABlineal[0];
        var b = this.getVariablesABlineal[1];
        var array = new Array();
        array = [a, b, c];
        return array;
    }
}

class LimpiarDatos {
    constructor(dat1, dat2) {
        this.__datX = dat1.split(',');
        this.__datY = dat2.split(',');
    }

    get getDataXdec() {
        let datosEnX = new Array();
        for (let i = 0; i < this.__datX.length; i++) {
            datosEnX[i] = parseFloat(this.__datX[i]);
        }
        return datosEnX;
    }

    get getDataYdec() {
        let datosEnY = new Array();
        for (let i = 0; i < this.__datY.length; i++) {
            datosEnY[i] = parseFloat(this.__datY[i]);
        }
        return datosEnY;
    }

    get getDatosXpot() {
        let datosEnX = new Array();
        var condicion = false;
        let menor = 1000000;

        for (let i = 0; i < this.__datX.length; i++) {
            if (parseFloat(this.__datX[i]) <= 0) {
                condicion = true;
                break;
            } else {
                continue;
            }
        }

        if (condicion == true) {
            for (let i = 0; i < this.__datX.length; i++) {
                if (parseFloat(this.__datX[i]) < menor) {
                    menor = parseFloat(this.__datX[i]);
                }
            }

            menor *= -1;

            for (let i = 0; i < this.__datX.length; i++) {
                datosEnX[i] = Math.log(parseFloat(this.__datX[i]) + menor + 1);
            }

            return datosEnX;
        } else {
            for (let i = 0; i < this.__datX.length; i++) {
                datosEnX[i] = Math.log(parseFloat(this.__datX[i]));
            }
            return datosEnX;
        }
    }

    get getDatosYpot() {
        let datosEnY = new Array();
        var condicion = false;
        let menor = 1000000;

        for (let i = 0; i < this.__datY.length; i++) {
            if (parseFloat(this.__datY[i]) <= 0) {
                condicion = true;
                break;
            } else {
                continue;
            }
        }

        if (condicion == true) {
            for (let i = 0; i < this.__datY.length; i++) {
                if (parseFloat(this.__datY[i]) < menor) {
                    menor = parseFloat(this.__datY[i]);
                }
            }

            menor *= -1;

            for (let i = 0; i < this.__datY.length; i++) {
                datosEnY[i] = Math.log(parseFloat(this.__datY[i]) + menor + 1);
            }

            return datosEnY;
        } else {
            for (let i = 0; i < this.__datY.length; i++) {
                datosEnY[i] = Math.log(parseFloat(this.__datY[i]));
            }
            return datosEnY;
        }
    }
}

var polin2 = new Array();
var funExp = new Array();
var funLin = new Array();
var funLog = new Array();
var funPot = new Array();

async function crearGrafica(x, y, numDatos) {
    const cleanData = new LimpiarDatos(x, y);

    const datosX = await cleanData.getDataXdec;
    const datosY = await cleanData.getDataYdec;
    const datosPotX = await cleanData.getDatosXpot;
    const datosPotY = await cleanData.getDatosYpot;

    if (datosX.length == datosY.length && datosX.length == numDatos) {
        document.getElementById('ocultar').style.display = "block";
        console.log('Datos X:');
        console.log(datosX);
        console.log('Datos Y:');
        console.log(datosY);
        crearTablaDeDatos(datosX, datosY, numDatos);
        var metodo = new MetodoPolinomial(datosX, datosY, numDatos);
        var potencia = new MetodoPolinomial(datosX, datosPotY, numDatos);
        var logaritmo = new MetodoPolinomial(datosPotX, datosY, numDatos);
        var ultimaFun = new MetodoPolinomial(datosPotX, datosPotY, numDatos);

        polin2 = metodo.getPolinomioResultanteGrado2;
        funExp = potencia.getValExp;
        funLin = metodo.getVariablesABlineal;
        funLog = logaritmo.getVariablesABlineal;
        funPot = ultimaFun.getVariablesABlineal;

        console.log(datosPotY);
        console.log(funExp);
        crearTablaMinimos(metodo.imprimirDatos, metodo.getSumasRows);
        addPolinomio(metodo.getPolinomioResultanteGrado2, metodo.getCorrelacion);
        graficaConPolinomios(datosX, datosY, metodo.getPolinomioResultanteGrado2, metodo.getVariablesABlineal, funExp, funLin, funLog, funPot);
        addLineal(metodo.getVariablesABlineal);
        return new Chart(document.getElementById("line-chart"), {
            type: 'line',
            data: {
                labels: datosX,
                datasets: [{
                    data: datosY,
                    label: "F(x)",
                    borderColor: "#3e95cd",
                    fill: false
                }]
            },
            options: {
                title: {
                    display: true,
                    text: 'Plano X - Y'
                }
            }
        });

    } else {
        alert(`debe ingresar ${numDatos} en las dos casillas datos en x: ${datosX.length}, datos en y: ${datosY.length}`);
        console.log(numDatos);
        console.log(datosY.length);
        console.log(datosX.length);
    }
}

function crearTablaDeDatos(x, y, numDatos) {
    const rowX = document.getElementById('rowChartX');

    var c1 = document.createElement('td');
    var contenido1 = document.createTextNode(`X`);
    c1.appendChild(contenido1);
    rowX.appendChild(c1);

    for (let index = 0; index < numDatos; index++) {
        var c = document.createElement('td');
        var contenido = document.createTextNode(x[index]);
        c.appendChild(contenido);
        rowX.appendChild(c);
    }

    const rowY = document.getElementById('rowChartY');

    var c2 = document.createElement('td');
    var contenido2 = document.createTextNode(`Y`);
    c2.appendChild(contenido2);
    rowY.appendChild(c2);

    for (let index = 0; index < numDatos; index++) {
        var c = document.createElement('td');
        var contenido = document.createTextNode(y[index]);
        c.appendChild(contenido);
        rowY.appendChild(c);
    }

}

function crearTablaMinimos(x, y) {
    const tabla = document.getElementById('tablaMinimos');

    for (let q = 0; q < x[0].length; q++) {
        var d1 = document.createElement('tr')
        for (let i = 0; i < 7; i++) {
            var c1 = document.createElement('td');
            var contenido1 = document.createTextNode(x[i][q]);
            c1.appendChild(contenido1);
            d1.appendChild(c1);
        }
        tabla.appendChild(d1);
    }

    const tabla2 = document.getElementById('tablaSumas');

    for (let i = 0; i < 7; i++) {
        var c = document.createElement('td');
        var con = document.createTextNode('Σ: ' + y[i]);
        c.appendChild(con);
        tabla2.appendChild(c);
    }

}

function addPolinomio(polinomio, corr) {
    var x1 = `${polinomio[2].toFixed(4)}x^2`;
    var x2 = '';
    var x3 = '';
    if (polinomio[1] >= 0) {
        x2 = `+${polinomio[1].toFixed(4)}x`;
    } else {
        x2 = `${polinomio[1].toFixed(4)}x`;
    }
    if (polinomio[0] >= 0) {
        x3 = `+${polinomio[0].toFixed(4)}`;
    } else {
        x3 = `${polinomio[0].toFixed(4)}`;
    }
    var c = document.createElement('h6');
    var con = document.createTextNode(`y(x)=${x1+x2+x3}`);
    c.appendChild(con);
    document.getElementById('contPolinomio').appendChild(c);

    var s = document.createElement('h6');
    var w = document.createTextNode(`r=${corr}`);
    s.appendChild(w);
    document.getElementById('contPolinomioR').appendChild(s);
}


async function graficaConPolinomios(xz1, yz1, pol1, ab, t, uu, u) {
    var mayor = 0;
    var menor = 10000000;
    var y = new Array();
    var x = new Array();
    var lineal = new Array();
    var expone = new Array();
    var loga = new Array();
    var pote = new Array();

    for (let index = 0; index < xz1.length; index++) {
        if (menor > xz1[index]) {
            menor = await xz1[index];
        }
    }

    for (let index = 0; index < xz1.length; index++) {
        if (mayor < xz1[index]) {
            mayor = await xz1[index];
        }
    }

    var cont = 0;
    var variable = menor;

    while (cont <= xz1.length) {
        expone[cont] = t[2] * Math.exp(t[0] * variable);
        loga[cont] = uu[0] + uu[1] * variable;
        pote[cont] = Math.exp(u[1]) * variable ^ u[0];
        y[cont] = await pol1[2] * xz1[cont] * xz1[cont];
        y[cont] += await pol1[1] * xz1[cont];
        y[cont] += await pol1[0];
        lineal[cont] = ab[0] * xz1[cont] + ab[1];
        x[cont] = await variable;
        cont++;
        variable += await (mayor - menor) / xz1.length;
    }

    return new Chart(document.getElementById("line-chart2"), {
        type: 'line',
        data: {
            labels: xz1,
            datasets: [{
                data: yz1,
                label: "puntos x - y",
                borderColor: "green",
                fill: false,
                yAxisID: 'y-axis-1'
            }, {
                data: y,
                label: "polinomio grado 2",
                borderColor: "red",
                fill: false,
                yAxisID: 'y-axis-2'
            }, {
                data: lineal,
                label: "Funcion Lineal",
                borderColor: "blue",
                fill: false,
                yAxisID: 'y-axis-3'
            }, {
                data: expone,
                label: "Funcion Exponencial",
                borderColor: "yellow",
                fill: false,
                yAxisID: 'y-axis-4'
            }, {
                data: loga,
                label: "Funcion Logaritmica",
                borderColor: "purple",
                fill: false,
                yAxisID: 'y-axis-5'
            }, {
                data: pote,
                label: "Funcion Potencial",
                borderColor: "grey",
                fill: false,
                yAxisID: 'y-axis-6'
            }]
        },
        options: {
            responsive: true,
            hoverMode: 'index',
            stacked: false,
            title: {
                display: true,
                text: 'graficas'
            },
            scales: {
                yAxes: [{
                    type: 'linear', // only linear but allow scale type registration. This allows extensions to exist solely for log scale for instance
                    display: true,
                    position: 'left',
                    id: 'y-axis-1',
                }, {
                    type: 'linear', // only linear but allow scale type registration. This allows extensions to exist solely for log scale for instance
                    display: true,
                    position: 'left',
                    id: 'y-axis-2',
                }, {
                    type: 'linear', // only linear but allow scale type registration. This allows extensions to exist solely for log scale for instance
                    display: true,
                    position: 'left',
                    id: 'y-axis-3',
                }, {
                    type: 'linear', // only linear but allow scale type registration. This allows extensions to exist solely for log scale for instance
                    display: true,
                    position: 'left',
                    id: 'y-axis-4',
                }, {
                    type: 'linear', // only linear but allow scale type registration. This allows extensions to exist solely for log scale for instance
                    display: true,
                    position: 'left',
                    id: 'y-axis-5',
                }, {
                    type: 'linear', // only linear but allow scale type registration. This allows extensions to exist solely for log scale for instance
                    display: true,
                    position: 'left',
                    id: 'y-axis-6',

                    // grid line settings
                    gridLines: {
                        drawOnChartArea: false, // only want the grid lines for one axis to show up
                    },
                }],
            }
        }
    });
}

function addLineal(ab) {
    if (ab[1] >= 0) {
        var lineal = ab[0].toFixed(4) + 'X' + `+${ab[1]}`;
    }
    var c = document.createElement('h6');
    var con = document.createTextNode(`y(x)=${lineal}`);
    c.appendChild(con);
    document.getElementById('contLineal').appendChild(c);
}

function calcInterpol(val) {
    var valor1 = polin2[2] * val * val + polin2[1] * val + polin2[0];
    var valor2 = funLin[0] * val + funLin[1];
    var valor3 = funExp[0] * Math.exp(funExp[1] * val);
    var valor4 = funLog[0] * Math.log(val) + funLog[1];
    var valor5 = Math.exp(funPot[1]) * val ^ funPot[0];
    var co = Math.exp(funPot[1]);
    var c = document.createElement('h6');
    var con = document.createTextNode(`funcion lineal: y(${val})=${valor2}`);
    c.appendChild(con);
    document.getElementById('interpolacionesCalculadas').appendChild(c);

    var x = document.createElement('h6');
    var qu = document.createTextNode(`polinomio grado 2: y(${val})=${valor1}`);
    x.appendChild(qu);
    document.getElementById('interpolacionesCalculadas').appendChild(x);

    var k = document.createElement('h6');
    var as = document.createTextNode(`Curva Exponencial: y(${val})=${funExp[2].toFixed(4)}e^(${funExp[0].toFixed(4)}x) = ${valor3.toFixed(4)}`);
    k.appendChild(as);
    document.getElementById('interpolacionesCalculadas').appendChild(k);

    var f = document.createElement('h6');
    var wwe = document.createTextNode(`Curva Logaritmica: y(${val})=${funLog[0].toFixed(4)}ln(x)+${funLog[1].toFixed(4)}=${valor4}`);
    f.appendChild(wwe);
    document.getElementById('interpolacionesCalculadas').appendChild(f);

    var rr = document.createElement('h6');
    var ert = document.createTextNode(`Curva Potencial: y(${val})=${co.toFixed(4)}*x^${funPot[0].toFixed(4)}=${valor5}`);
    rr.appendChild(ert);
    document.getElementById('interpolacionesCalculadas').appendChild(rr);

    var rr = document.createElement('hr');
    document.getElementById('interpolacionesCalculadas').appendChild(rr);
}

/*
console.log('Ecuaciones Sistema 3*3');
console.table(metodo1.getPolinomio);
console.log('Sistema 3*3');
console.table(metodo1.getMatriz);
console.log('Objetivos');
console.log('Polinomio Grado 2: ');
console.table(`${metodo1.getPolinomioResultanteGrado2[2].toFixed(4)}X^2 + ${metodo1.getPolinomioResultanteGrado2[1].toFixed(4)}X + ${metodo1.getPolinomioResultanteGrado2[0].toFixed(4)}`);
console.log('polinomio interpolado en x = 0.7');
console.log(`y(0.7) = ${metodo1.getImagenPolinomioGrado2}`)
-1,0,1,2,3
6.62,3.94,2.17,1.35,0.89
*/