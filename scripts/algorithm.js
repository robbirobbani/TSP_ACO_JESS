function ACO(data,t1) {
  var runningtime2 = new Array(100)
  runningtime2[0] = 0
  var Rho = 0.5;
  var Thao = 0.2;
  var alpha = 1;
  var betha = 1;
  var m = 5;
  var n = 5;
  var lenD = data.length;
  var Iter = 0;

  var maxIter = 100;

  var temp = new Array(data.length + 1);
  var posisi = randPosisi(data, 1);

  // console.log(data);

  var randomSemut = pengkodeanPopulasi(posisi);
  // console.log(randomSemut);
  var urutandummy = [[1, 0, 3, 4, 2]];
  // console.log(urutandummy);

  var penyebut = nilaiPenyebut(Thao, alpha, betha, data);
   console.log(penyebut);
 

  var kelilingKota = urutanKota(
    Thao,
    alpha,
    betha,
    penyebut,
    data,
    randomSemut
  );
   console.log(kelilingKota);

  var panjang = panjangPerjalanan(kelilingKota, data);
  // console.log(panjang);

  var ThaoBaru = newThao(Rho, Thao);
  // console.log(ThaoBaru);

  var matrixThao = matrixNewT(kelilingKota, panjang, ThaoBaru);
  // console.log(matrixThao);

  var matrixThaoSim = matrixTSim(matrixThao);
  // console.log(matrixThaoSim);
  // console.log(kelilingKota);
  var savePanjang = saveArray(kelilingKota, Iter, temp, panjang);
  // console.log(savePanjang);
  runningtime2[1] = performance.now()-t1
  Iter = 1;

  var penyebut = nilaiPenyebut(Thao, alpha, betha, data);
  // console.log(penyebut);
  

  var kelilingKota = urutanKota(
    Thao,
    alpha,
    betha,
    penyebut,
    data,
    randomSemut
  );
  // console.log(kelilingKota);

  var panjang = panjangPerjalanan(kelilingKota, data);
  // console.log(panjang);

  var ThaoBaru = newThao(Rho, Thao);
  // console.log(ThaoBaru);

  var matrixThao = matrixNewT(kelilingKota, panjang, ThaoBaru);
  // console.log(matrixThao);

  var matrixThaoSim = matrixTSim(matrixThao);
  // console.log(matrixThaoSim);

  var newIterThao = newThao2(matrixThaoSim, Rho);
  // console.log(newIterThao);

  savePanjang = saveArray(kelilingKota, Iter, temp, panjang);
  // console.log(savePanjang);
  runningtime2[2] = performance.now()-t1;
  for (var i = 2; i < maxIter; i++) {
    Iter = i;
    var posisi = randPosisi(data, 1);

    // console.log(data);

    var randomSemut = pengkodeanPopulasi(posisi);
    // console.log(randomSemut);

    var penyebut = nilaiPenyebut(Thao, alpha, betha, data);
    // console.log(penyebut);
   

    var kelilingKota = urutanKota(
      Thao,
      alpha,
      betha,
      penyebut,
      data,
      randomSemut
    );
    // console.log(kelilingKota);

    var panjang = panjangPerjalanan(kelilingKota, data);
    // console.log(panjang);

    // var ThaoBaru = newThao(Rho, Thao);
    // console.log(ThaoBaru);

    var matrixThao = matrixNewT2(kelilingKota, panjang, newIterThao);
    // console.log(matrixThao);

    var matrixThaoSim = matrixTSim(matrixThao);
    // console.log(matrixThaoSim);

    var newIterThao = newThao2(matrixThaoSim, Rho);
    // console.log(newIterThao);

    savePanjang = saveArray(kelilingKota, Iter, temp, panjang);
    // console.log(savePanjang);
    console.log(Iter);
    runningtime2[Iter+1]= performance.now()-t1
  }

  var hasil = hasilTSP(data, maxIter, savePanjang);
  // console.log(hasil);

  var hasilJarak = hasil[data.length];
  console.log(hasilJarak);
  var hasilIterasi = hasil[data.length + 1] + 1;
  console.log(hasilIterasi);

  hasil.pop();
  hasil.pop();
  // console.log(hasil);

  for (var i = 0; i < data.length; i++) {
    hasil[i] += 1;
  }

  console.log(hasil);
  console.log(runningtime2)
  // var Random =
  // var fitness = hitungFitnessPopulasi(data,randomSemut);						//Fitness
  // console.log(fitness);
}

function MACO(data,t0) {
  var runningtime =  new Array(100)
  runningtime [0] =  0
 
  var Rho = 0.5;
  var Thao = 0.2;
  var alpha = 1;
  var betha = 1;
  var m = 5;
  var n = 5;
  var lenD = data.length;
  var Iter = 0;

  var maxIter = 100;

  var temp = new Array(data.length + 1);
  var posisi = randPosisi(data, 1);

  // console.log(data);

  var randomSemut = pengkodeanPopulasi(posisi);

  var bestTemp = thaoMACO(data, randomSemut);
  // console.log(randomSemut);
  var urutandummy = [[1, 0, 3, 4, 2]];
  // console.log(urutandummy);

  var penyebut = nilaiPenyebut(Thao, alpha, betha, data);
  // console.log(penyebut);
  var sortRS = sortPenyebut(randomSemut, penyebut);
  // console.log(sortRS);

  var kelilingKota = urutanKotaMACO(
    Thao,
    alpha,
    betha,
    penyebut,
    data,
    randomSemut,
    bestTemp
  );
  // console.log(kelilingKota);

  var panjang = panjangPerjalanan(kelilingKota, data);
  // console.log(panjang);

  var ThaoBaru = newThao(Rho, Thao);
  // console.log(ThaoBaru);

  var matrixThao = matrixNewT(kelilingKota, panjang, ThaoBaru);
  // console.log(matrixThao);

  var matrixThaoSim = matrixTSim(matrixThao);
  // console.log(matrixThaoSim);
  // console.log(kelilingKota);
  var savePanjang = saveArray(kelilingKota, Iter, temp, panjang);
  // console.log(savePanjang);
  runningtime[1]=performance.now()-t0
  Iter = 1;

  var penyebut = nilaiPenyebut(Thao, alpha, betha, data);
  // console.log(penyebut);
  var sortRS = sortPenyebut(randomSemut, penyebut);
  // console.log(sortRS);

  var kelilingKota = urutanKotaMACO(
    Thao,
    alpha,
    betha,
    penyebut,
    data,
    randomSemut,
    bestTemp
  );
  // console.log(kelilingKota);

  var panjang = panjangPerjalanan(kelilingKota, data);
  // console.log(panjang);

  var ThaoBaru = newThao(Rho, Thao);
  // console.log(ThaoBaru);

  var matrixThao = matrixNewT(kelilingKota, panjang, ThaoBaru);
  // console.log(matrixThao);

  var matrixThaoSim = matrixTSim(matrixThao);
  // console.log(matrixThaoSim);

  var newIterThao = newThao2(matrixThaoSim, Rho);
  // console.log(newIterThao);

  savePanjang = saveArray(kelilingKota, Iter, temp, panjang);
  // console.log(savePanjang);
  runningtime[2]= performance.now()-t0

  for (var i = 2; i < maxIter; i++) {
    Iter = i;
    var posisi = randPosisi(data, 1);

    // console.log(data);

    var randomSemut = pengkodeanPopulasi(posisi);

    // console.log(randomSemut);

    var penyebut = nilaiPenyebut(Thao, alpha, betha, data);
    // console.log(penyebut);
    var sortRS = sortPenyebut(randomSemut, penyebut);
    // console.log(sortRS);

    var kelilingKota = urutanKotaMACO(
      Thao,
      alpha,
      betha,
      penyebut,
      data,
      randomSemut,
      bestTemp
    );
    // console.log(kelilingKota);

    var panjang = panjangPerjalanan(kelilingKota, data);
    // console.log(panjang);

    // var ThaoBaru = newThao(Rho, Thao);
    // console.log(ThaoBaru);

    var matrixThao = matrixNewT2(kelilingKota, panjang, newIterThao);
    // console.log(matrixThao);

    var matrixThaoSim = matrixTSim(matrixThao);
    // console.log(matrixThaoSim);

    var newIterThao = newThao2(matrixThaoSim, Rho);
    // console.log(newIterThao);

    savePanjang = saveArray(kelilingKota, Iter, temp, panjang);
    // console.log(savePanjang);
    runningtime[Iter+1]= performance.now()-t0
    console.log(Iter);
  }

  var hasil = hasilTSP(data, maxIter, savePanjang);
  // console.log(hasil);

  var hasilJarak = hasil[data.length];
  console.log(hasilJarak);
  var hasilIterasi = hasil[data.length + 1] + 1;
  console.log(hasilIterasi);

  
  // console.log(hasil);

  for (var i = 0; i < data.length; i++) {
    hasil[i] += 1;
  }

  console.log(hasil);
  console.log(bestTemp);
  console.log(runningtime)
  // var Random =
  // var fitness = hitungFitnessPopulasi(data,randomSemut);           //Fitness
  // console.log(fitness);
}

//Copy Matrix 2-dimension
function copyD(d) {
  var copyData = new Array();
  for (var i = 0; i < d.length; i++) {
    copyData[i] = [];
    for (var j = 0; j < d[0].length; j++) {
      copyData[i][j] = d[i][j];
    }
  }
  return copyData;
}
//Copy Vector
function copyV(v) {
  var copyVec = new Array();
  for (var i = 0; i < v.length; i++) {
    copyVec[i] = v[i];
  }
  return copyVec;
}
//Copy Var
function copyVar(v) {
  var copyVar = v;
  return v;
}
//Random Posisi Partikel Awal
function randPosisi(data, par) {
  var randomPos = new Array();
  for (var i = 0; i < par; i++) {
    randomPos[i] = [];
    for (var j = 0; j < data.length; j++) {
      randomPos[i][j] = Math.random();
    }
  }
  return randomPos;
}
//Random Kecepatan Awal
function randKecepatan(data, par) {
  var randomKec = new Array();
  for (var i = 0; i < par; i++) {
    randomKec[i] = [];
    for (var j = 0; j < data.length; j++) {
      randomKec[i][j] = Math.random();
    }
  }
  return randomKec;
}
//Konversi pengkodean TSP populasi
function pengkodeanPopulasi(pos) {
  var urutan = new Array();
  //Pencocokan Posisi
  for (var i = 0; i < pos.length; i++) {
    urutan[i] = pengkodean1Partikel(pos[i]);
  }
  return urutan;
}

function pengkodean1Partikel(par) {
  var dummy = copyV(par);
  var urutan = new Array();
  //Sorting
  dummy.sort();
  //Pencocokan Posisi
  for (var i = 0; i < par.length; i++) {
    for (var j = 0; j < par.length; j++) {
      if (par[i] == dummy[j]) {
        urutan[i] = j;
      }
    }
  }
  return urutan;
}

function nilaiPenyebut(th, a, b, data) {
  var np = new Array();
  var sum = 0;
  var x = 0;
  var y = 0;
  var z = 0;
  var tot = 0;

  for (var i = 0; i < data.length; i++) {
    np[i] = [];
    for (var j = 0; j < data.length; j++) {
      if (i != j) {
        // x = Math.pow(th,a);
        // y = Math.pow((1/data[i][j]),b)
        sum = Math.pow(th, a) * Math.pow(1 / data[i][j], b);
        tot += sum;
        // console.log(tot)
      }
    }
    np[i][1] = tot;
    np[i][0] = i;
    tot = 0;
  }

  return np;
}

function sortPenyebut(rs, np) {
  var copyRS = new Array();
  for (var i = 0; i < np.length; i++) {
    copyRS[i] = [];
    for (var j = 0; j < np.length; j++) {
      if (rs[0][i] == np[j][0]) {
        copyRS[i][1] = np[j][1];
        copyRS[i][0] = rs[0][i];
      }
    }
  }
  return copyRS;
}

function urutanKota(th, a, b, copyRS, data, randurt) {
  // console.log(copyRS);
  var jalur = [];
  var jalurKota = [];
  var bilrandom = [];
  var sum = 0;
  var sum2 = 0;
  var sum3 = 0;
  var sum4 = 0;
  var cek = 0;

  for (var i = 0; i < data.length; i++) {
    jalur[i] = [];

    bilrandom[i] = [];
    for (var j = 0; j < data.length; j++) {
      if (j == randurt[0][i]) {
        sum = 0;
        jalur[i][j] = sum;
      } else {
        sum =
          (Math.pow(th, a) * Math.pow(1 / data[randurt[0][i]][j], b)) /
          copyRS[randurt[0][i]][1];

        jalur[i][j] = sum;
      }
    }
  }

  for (var i = 0; i < data.length; i++) {
    jalurKota[i] = [];
    jalurKota[data.length] = [];
    jalurKota[0][i] = randurt[0][i];
  }

  for (var h = 1; h < data.length; h++) {
    // console.log(jalurKota);

    for (var i = 0; i < data.length; i++) {
      for (var m = 0; m < data.length; m++) {
        if (jalur[i][m] != 0) {
          var tmprandom = Math.random();
          bilrandom[i] = tmprandom;
          // console.log(jalur[i][m]);
          // console.log(bilrandom[i]);

          if (jalur[i][m] > bilrandom[i]) {
            bilrandom[i] = bilrandom[i] + jalur[i][m];
            // console.log(bilrandom[i]);
            m = m + data.length;
          }
          m = m + data.length;
        }
      }
      // console.log(bilrandom);

      for (var j = 0; j < data.length; j++) {
        // console.log("AA" + i + "dan" + j);
        if (
          jalur[i][j] != 0 &&
          sum2 + jalur[i][j] < bilrandom[i] &&
          j != data.length - 1
        ) {
          cek = j;
          sum2 = sum2 + jalur[i][j];
          // console.log("SATU");
          // console.log(jalur[i][j]);
        } else if (sum2 + jalur[i][j] > bilrandom[i] && j != data.length - 1) {
          // console.log(jalur[i][cek]);
          // console.log(cek);
          jalurKota[h][i] = cek;
          sum2 = 0;
          jalur[i][cek] = sum2;
          // console.log(jalur[i][cek]);
          j = j + data.length;
          // console.log("DUA");
          // console.log(jalur);
        } else if (j == data.length - 1 && jalur[i][j] != 0) {
          cek = j;
          jalur[i][j] = 0;
          jalurKota[h][i] = cek;
          sum2 = 0;
          // console.log("TIGA");
          // console.log(jalur[i][j]);
        } else if (j == data.length - 1 && jalur[i][j] == 0) {
          jalur[i][cek] = 0;
          jalurKota[h][i] = cek;
          sum2 = 0;
          // console.log("EMPAT");
          // console.log(jalur[i][cek]);
        } else {
        }
        // console.log("LIMA");
      }

      for (var k = 0; k < data.length; k++) {
        if (jalur[i][k] != 0) {
          // console.log(jalur[i][k] + "HH" + k);
          // console.log(copyRS[jalurKota[h][i]][1]);
          // console.log([jalurKota[h][i]] + "dan" + k + "dan" + h);
          // console.log(data[jalurKota[h][i]][k]);
          sum3 =
            (Math.pow(th, a) * Math.pow(1 / data[jalurKota[h][i]][k], b)) /
            copyRS[jalurKota[h][i]][1];
          jalur[i][k] = sum3;
          // console.log(jalur[i][k] + "JJ" + k);
        }
      }
    }
  }

  // console.log(bilrandom);
  // console.log(jalurKota);

  // console.log(jalur);
  for (var i = 0; i < data.length; i++) {
    jalurKota[data.length][i] = jalurKota[0][i];
  }
  return jalurKota;
}

function panjangPerjalanan(jalurkota, data) {
  var total = new Array();
  var sum = 0;
  // console.log(jalurkota);
  for (var j = 0; j < data.length; j++) {
    total[j] = [];
    for (var k = 0; k < data.length; k++) {
      // console.log(data[jalurkota[k][j]]);
      // console.log(jalurkota[k + 1][j]);
      sum += data[jalurkota[k][j]][jalurkota[k + 1][j]];
      total[j] = sum;
    }
    sum = 0;
  }
  return total;
}

function newThao(rho, thao) {
  var NT = 0;
  NT = (1 - rho) * thao;

  return NT;
}

function matrixNewT(jalurkota, panjang, NT) {
   var newThao = new Array(jalurkota.length - 1);
  var sum = 0;
  
    for (var j = 0; j < jalurkota.length - 1; j++) {
      newThao[j]=[];
      newThao[j][j] = 0;
    }
 for (var l = 0; l < jalurkota.length - 1; l++) {
        for (var k = 0; k < jalurkota.length - 1; k++) {
          if (k != l && jalurkota[k][l] == k && jalurkota[k + 1][l] == l) {
           
            sum += 1 / panjang[l];
          
            newThao[k][l] = sum;
          }
        }
        sum = 0;
      }
  for (var i = 0; i < jalurkota.length - 1; i++) {
    for (var j = 0; j < jalurkota.length - 1; j++) {
      if (i != j) {
        newThao[i][j] += NT;
      }
    }
  }
  return newThao;
}

function matrixTSim(NT) {
  var simThao = new Array();
  for (var i = 0; i < NT.length; i++) {
    simThao[i] = [];
    for (var j = 0; j < NT.length; j++) {
      simThao[i][j] = 0;
    }
  }
  for (var i = 0; i < NT.length; i++) {
    for (var j = 0; j < NT.length; j++) {
      if (i < j) {
        simThao[i][j] = NT[i][j] + NT[j][i];
        simThao[j][i] = simThao[i][j];
      }
    }
  }
  return simThao;
}

function newThao2(simThao, Rho) {
  var newNT = new Array();
  for (var i = 0; i < simThao.length; i++) {
    newNT[i] = [];
    for (var j = 0; j < simThao.length; j++) {
      newNT[i][j] = 0;
    }
  }
  for (var k = 0; k < simThao.length; k++) {
    for (var l = 0; l < simThao.length; l++) {
      if (k != l) newNT[k][l] = (1 - Rho) * simThao[k][l];
    }
  }
  return newNT;
}

function matrixNewT2(jalurkota, panjang, NT) {
  var newThao = new Array(jalurkota.length - 1);
  var sum = 0;
  
    
    for (var j = 0; j < jalurkota.length - 1; j++) {
      newThao[j] = [];
      newThao[j][j] = 0;
    }

  for (var l = 0; l < jalurkota.length - 1; l++) {
        for (var k = 0; k < jalurkota.length - 1; k++) {
          if (k != l && jalurkota[k][l] == k && jalurkota[k + 1][l] == l) {
            sum += 1 / panjang[l];
            newThao[k][l] = sum;
          }
        }
        sum = 0;
      }

  for (var i = 0; i < jalurkota.length - 1; i++) {
    for (var j = 0; j < jalurkota.length - 1; j++) {
      if (i != j) {
        newThao[i][j] += NT[i][j];
      }
    }
  }
  return newThao;
}

function saveArray(d, Iter, temp, panjang) {
  temp[Iter] = [];
  for (var i = 0; i < d.length; i++) {
    temp[Iter][i] = [];
    for (var j = 0; j < d.length - 1; j++) {
      temp[Iter][i][j] = d[i][j];

      // copyData[i][j] = d[i][j];
    }
  }
  for (var i = 0; i < d.length - 1; i++) {
    temp[Iter][d.length - 1][i] = panjang[i];
  }
  return temp;
}

function hasilTSP(data, maxIter, savePanjang) {
  var hasil = new Array(data.length + 2);
  for (var i = 0; i < data.length + 2; i++) {
    hasil[i] = [];
  }
  for (var i = 0; i < maxIter; i++) {
    // savePanjang[i]
    // hasil[i]  = [];
    for (var j = 0; j < data.length; j++) {
      // savePanjang[i]
      if (i == 0 && j == 0) {
        hasil[data.length] = savePanjang[0][data.length][0];
        hasil[data.length + 1] = i;
        for (var k = 0; k < data.length; k++) {
          hasil[k] = savePanjang[i][k][j];
          // console.log(hasil)
        }
      } else if (savePanjang[i][data.length][j] < hasil[data.length]) {
        hasil[data.length] = savePanjang[i][data.length][j];
        hasil[data.length + 1] = i;
        for (var k = 0; k < data.length; k++) {
          hasil[k] = savePanjang[i][k][j];
        }
      }
      // hasil[i][j]
    }
  }
  return hasil;
}


  function thaoMACO(data, urutan){
  var dataThao = new Array();
  var dataX = new Array();
  var Thao = new Array();
  dataX=copyD(data)
  console.log(dataX)
  for (var i = 0; i < data.length; i++) {
    
    dataX[i][i] = 9007199254740992;
  } console.log(Thao)
  dataThao = copyD(dataX);
  console.log(dataThao)
  for (var i = 0; i < data.length; i++) {
    Thao[i] = [];
    Thao[i][data.length] = [];

       for (var j = 0; j < data.length; j++) {
      if(j==0){
        
        Thao[i][j] = urutan[0][i]
        
        var min = Math.min.apply(Math, dataThao[Thao[i][j]]);
        
        Thao[i][j+1] = dataThao[Thao[i][j]].indexOf(min)
        
        Thao[i][data.length] = min
        for (var k = 0; k < data.length; k++) {
          dataThao[k][urutan[0][i]] = 9007199254740992;
        }
        for (var k = 0; k < data.length; k++) {
          dataThao[k][Thao[i][j+1]] = 9007199254740992;
        }
      } else  if(j!=data.length-1){
        
        var min = Math.min.apply(Math, dataThao[Thao[i][j]]);
       
        Thao[i][j+1] = dataThao[Thao[i][j]].indexOf(min)
        
        Thao[i][data.length] += min
        for (var l = 0; l < data.length; l++) {

          dataThao[l][Thao[i][j]] = 9007199254740992;
        }
        
      } else {
        var min = dataX[Thao[i][j]][Thao[i][0]]
        Thao[i][data.length] += min
      }

    }
    dataThao= copyD(dataX);

  }
  
  var jalurterkecil = new Array(data.length+3)
  jalurterkecil[data.length+1]= Thao[0][data.length]
  jalurterkecil[data.length+2]= 0
  for (var i = 0; i < data.length; i++) {
    if(jalurterkecil[data.length+1]>Thao[i][data.length]){
      
        jalurterkecil[data.length+1]= Thao[i][data.length]
        jalurterkecil[data.length+2]= i
    }
  }
  console.log(jalurterkecil)
  for (var i = 0; i < data.length; i++) {
    
        jalurterkecil[i]= Thao[jalurterkecil[data.length+2]][i]
  }
  jalurterkecil[data.length]=jalurterkecil[0];
  jalurterkecil[data.length+2]=1/(data.length*jalurterkecil[data.length+1]);
  


return jalurterkecil;
}

function urutanKotaMACO(th, a, b, copyRS, data, randurt,jalurterkecil) {
  // console.log(copyRS);
  var jalur = [];
  var jalurKota = [];
  var bilrandom = [];
  var sum = 0;
  var sum2 = 0;
  var sum3 = 0;
  var sum4 = 0;
  var cek = 0;

  for (var i = 0; i < data.length; i++) {
    jalur[i] = [];

    bilrandom[i] = [];
    for (var j = 0; j < data.length; j++) {
      if (j == randurt[0][i]) {
        sum = 0;
        jalur[i][j] = sum;
      } else {
        sum = Math.pow(1 / data[randurt[0][i]][j], b)
        if(sum>jalurterkecil[data.length+2]){       
          jalur[i][j] = sum; 
        }
        else  { jalurterkecil[i][j] = jalurterkecil[data.length+2]
        }
      }
    }
  }

  for (var i = 0; i < data.length; i++) {
    jalurKota[i] = [];
    jalurKota[data.length] = [];
    jalurKota[0][i] = randurt[0][i];
  }

  for (var h = 1; h < data.length; h++) {
   

    for (var i = 0; i < data.length; i++) {
      for (var m = 0; m < data.length; m++) {
        if (jalur[i][m] != 0) {
          var tmprandom = Math.random();
          bilrandom[i] = tmprandom;
          

          if (jalur[i][m] > bilrandom[i]) {
            bilrandom[i] = bilrandom[i] + jalur[i][m];
            
            m = m + data.length;
          }
          m = m + data.length;
        }
      }
   

      for (var j = 0; j < data.length; j++) {
        
        if (
          jalur[i][j] != 0 &&
          sum2 + jalur[i][j] < bilrandom[i] &&
          j != data.length - 1
        ) {
          cek = j;
          sum2 = sum2 + jalur[i][j];
         
        } else if (sum2 + jalur[i][j] > bilrandom[i] && j != data.length - 1) {
         
          jalurKota[h][i] = cek;
          sum2 = 0;
          jalur[i][cek] = sum2;
         
          j = j + data.length;
        
        } else if (j == data.length - 1 && jalur[i][j] != 0) {
          cek = j;
          jalur[i][j] = 0;
          jalurKota[h][i] = cek;
          sum2 = 0;
          
        } else if (j == data.length - 1 && jalur[i][j] == 0) {
          jalur[i][cek] = 0;
          jalurKota[h][i] = cek;
          sum2 = 0;
          
        } else {

        }
       
      }

      for (var k = 0; k < data.length; k++) {
        if (jalur[i][k] != 0) {
          
          sum = Math.pow(1 / data[randurt[0][i]][k], b)
        if(sum>jalurterkecil[data.length+2]){       
          jalur[i][k] = sum; 
        }
        else  { jalurterkecil[i][k] = jalurterkecil[data.length+2]
        }
         
        }
      }
    }
  }

 
  for (var i = 0; i < data.length; i++) {
    jalurKota[data.length][i] = jalurKota[0][i];
  }
  return jalurKota;
}
