// var wMax=0.9;
// var wMin=0.05;
// var c1=2;
// var c2=2;
// var r1;
// var r2;

function ACO(data) {
  var Rho = 0.5;
  var Thao = 0.2;
  var alpha = 1;
  var betha = 1;
  var m = 5;
  var n = 5;
  var lenD = data.length;
  var posisi = randPosisi(data, 1);

  console.log(data);

  var randomSemut = pengkodeanPopulasi(posisi);
  console.log(randomSemut);
  var urutandummy = [[1, 0, 3, 4, 2]];
  console.log(urutandummy);

  var penyebut = nilaiPenyebut(Thao, alpha, betha, data);
  console.log(penyebut);
  var sortRS = sortPenyebut(urutandummy, penyebut);
  console.log(sortRS);

  var kelilingKota = urutanKota(
    Thao,
    alpha,
    betha,
    penyebut,
    data,
    urutandummy
  );
  console.log(kelilingKota);

  // var Random =
  // var fitness = hitungFitnessPopulasi(data,randomSemut);						//Fitness
  // console.log(fitness);
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
  console.log(copyRS);
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

    jalurKota[0][i] = randurt[0][i];
  }

  for (var h = 1; h < data.length; h++) {
    console.log(jalurKota);

    for (var i = 0; i < data.length; i++) {
      for (var m = 0; m < data.length; m++) {
        if (jalur[i][m] != 0) {
          var tmprandom = Math.random();
          bilrandom[i] = tmprandom;
          console.log(jalur[i][m]);
          console.log(bilrandom[i]);

          if (jalur[i][m] > bilrandom[i]) {
            bilrandom[i] = bilrandom[i] + jalur[i][m];
            console.log(bilrandom[i]);
            m = m + data.length;
          }
          m = m + data.length;
        }
      }
      console.log(bilrandom);

      for (var j = 0; j < data.length; j++) {
        console.log("AA" + i + "dan" + j);
        if (
          jalur[i][j] != 0 &&
          sum2 + jalur[i][j] < bilrandom[i] &&
          j != data.length - 1
        ) {
          cek = j;
          sum2 = sum2 + jalur[i][j];
          console.log("SATU");
          console.log(jalur[i][j]);
        } else if (sum2 + jalur[i][j] > bilrandom[i] && j != data.length - 1) {
          console.log(jalur[i][cek]);
          console.log(cek);
          jalurKota[h][i] = cek;
          sum2 = 0;
          jalur[i][cek] = sum2;
          console.log(jalur[i][cek]);
          j = j + data.length;
          console.log("DUA");
          console.log(jalur);
        } else if (j == data.length - 1 && jalur[i][j] != 0) {
          cek = j;
          jalur[i][j] = 0;
          jalurKota[h][i] = cek;
          sum2 = 0;
          console.log("TIGA");
          console.log(jalur[i][j]);
        } else if (j == data.length - 1 && jalur[i][j] == 0) {
          jalur[i][cek] = 0;
          jalurKota[h][i] = cek;
          sum2 = 0;
          console.log("EMPAT");
          console.log(jalur[i][cek]);
        } else console.log("LIMA");
      }

      for (var k = 0; k < data.length; k++) {
        if (jalur[i][k] != 0) {
          console.log(jalur[i][k] + "HH" + k);
          console.log(copyRS[jalurKota[h][i]][1]);
          console.log([jalurKota[h][i]] + "dan" + k + "dan" + h);
          console.log(data[jalurKota[h][i]][k]);
          sum3 =
            (Math.pow(th, a) * Math.pow(1 / data[jalurKota[h][i]][k], b)) /
            copyRS[jalurKota[h][i]][1];
          jalur[i][k] = sum3;
          console.log(jalur[i][k] + "JJ" + k);
        }
      }
    }
  }

  console.log(bilrandom);
  console.log(jalurKota);

  console.log(jalur);
  return jalur;
}

// function acakSemut()
