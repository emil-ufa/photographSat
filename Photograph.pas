﻿unit Photograph;

interface
uses
    Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
    Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, VCLTee.TeEngine,
    VCLTee.Series, Vcl.ExtCtrls, VCLTee.TeeProcs, VCLTee.Chart, Vcl.ComCtrls, System.DateUtils,
    Bal_Types, anomal, astronom, deftorm, numint, atmos, CatDB, precpred, pg_type, math,
    NoradDB, orbint, sgp_h, Lagrange, VCLTee.GanttCh, prognozt, VclTee.TeeGDIPlus;

type
    TForm4 = class(TForm)
    Button1: TButton;
    Chart1: TChart;
    Series1: TPointSeries;
    maxDistance: TEdit;
    maxDistanceLabel: TLabel;
    satelliteNumberLabel: TLabel;
    satelliteNumber: TEdit;
    timeSampleLabel: TLabel;
    timeSample: TEdit;
    predictionTimeLabel: TLabel;
    predictionTime: TEdit;
    Chart2: TChart;
    PointSeries1: TPointSeries;
    clearDiagrams: TCheckBox;
    currentPredictionDay: TLabel;
        {Запускает программу построения диаграммы}
        procedure Button1Click(Sender: TObject);
        private
            { Private declarations }
        public
            { Public declarations }
    end;

    TSatOrb = record
        JDDouble: Double;
        Kdelta: Double;
        orbL: TLagrange;
        orbX: TXYZVxVyVz;
        case integer of
            0 : (orbK: TIOrb;);
            1 : (items : array [0..5] of Double;);
        end;

    function distanceBetweenSatellites(x1, x2: TXYZVxVyVz; dt : Double): Double;
    procedure minDistanceBetweenOrbits(orb1, orb2 : TIOrb; var r11, r12, r21, r22, p11, p12, p21, p22: Double);
    procedure readConfigVarsFromScreen(dt, predictionDays : integer; maxDist : Double; form: TForm4);
    function initSatDatabase(pathToOrbFiles, fileNameORB: string) : array of TSatOrb;
var
    Form4: TForm4;

implementation
{$R *.dfm}

// Функция рассчитывает расстояние между двумя КО, положение которых измеряется
// каждые dt секунд. Если за время dt достигалось минимальное расстояние между
// КО (расстояние от одного из них до прямой, по которой летит другой),
// тогда функция возвращает это минимальное расстояние.
//
// Входные параметры:
//      x1, x2 - векторы состояния каждого КО
//      dt     - временной интервал измерения положения КО
function distanceBetweenSatellites(x1, x2: TXYZVxVyVz; dt : Double): Double;
var
    r1, r2, v1, v2, dv, dri, drf : TXYZ;
    i                            : integer;
    abs_dv, abs_dri, abs_drf     : Double;
    dri_dv_scal, drf_dv_scal     : Double;
    cosi, cosf                   : Double;
begin
    r1 := x1.pos; v1 := x1.vel;
    r2 := x2.pos; v2 := x2.vel;

    dri_dv_scal := 0; drf_dv_scal := 0;

    for i := 0 to 2 do begin
        dri.items[i] := r2.items[i] - r1.items[i];
        dv.items[i] := v1.items[i] - v2.items[i];
        drf.items[i] := dri.items[i] - dv.items[i]*dt;
        dri_dv_scal := dri_dv_scal + dri.items[i] * dv.items[i];
        drf_dv_scal := drf_dv_scal + drf.items[i] * dv.items[i];
    end;

    abs_dri := sqrt(sqr(dri.x) + sqr(dri.y) + sqr(dri.z));
    abs_drf := sqrt(sqr(drf.x) + sqr(drf.y) + sqr(drf.z));
    abs_dv := sqrt(sqr(dv.x) + sqr(dv.y) + sqr(dv.z));

    cosi := dri_dv_scal/(abs_dri * abs_dv);
    cosf := drf_dv_scal/(abs_drf * abs_dv);

    if (cosi > 0) and (cosf <= 0) then begin
        result := sqrt(sqr(dri.y*dv.z - dri.z*dv.y) +
                       sqr(dri.z*dv.x - dri.x*dv.z) +
                       sqr(dri.x*dv.y - dri.y*dv.x))/abs_dv;
    end else begin
        result := abs_dri;
    end;
end;

// Функция рассчитывает минимальное расстояние между орбитами
// двух КО по линии узлов (ЛУ) - прямой пересечения плоскостей орбит.
//
// Входные параметры:
//      orb1, orb2 - векторы состояния двух КО (в элементах Кеплера)
//
// Выходные параметры:
//      r11, r12, r21, r22 - расстояния от притягивающего центра до точек пересечения орбит с ЛУ
//      p11, p12, p21, p22 - аргументы точек пересечения орбит с ЛУ
procedure minDistanceBetweenOrbits(orb1, orb2 : TIOrb; var r11, r12, r21, r22, p11, p12, p21, p22: Double);
var
    X : TXYZVxVyVz; // промежуточный вектор в декартовых координатах
    orbK: TIOrb; // промежуточные элементы Кеплера
    orbL : TLagrange; // промежуточные элементы Лагранжа
    ap_vector1, ap_vector2 : TXYZ; // радиус-векторы перигеев орбит
    common_vector : TXYZ; // направляющий вектор линии узлов
    second_ap_vector1, second_ap_vector2: TXYZ; // векторы, задающие соответственно плоскость 1й и 2й орбиты
    n1, n2: TXYZ; // нормали к плоскостям орбит
    cosf1, sinf1, cosf2, sinf2 : Double; // cos и sin углов, задающих точки пересечения орбит с линией узлов
    abs_common_vector, abs_ap_vector1, abs_ap_vector2 : Double; // модули используемых векторов
    phi : Double; // промежуточный угол
begin
    orbK := orb1;
    orbK.v := 0;
    keplerToGNSK(orbK, X);
    ap_vector1 := X.pos;

    orbK := orb1;
    orbK.v := Pi/2;
    keplerToGNSK(orbK, X);
    second_ap_vector1 := X.pos;

    orbK := orb2;
    orbK.v := 0;
    keplerToGNSK(orbK, X);
    ap_vector2 := X.pos;

    orbK := orb2;
    orbK.v := Pi/2;
    keplerToGNSK(orbK, X);
    second_ap_vector2 := X.pos;

    n1 := crossProduct(ap_vector1, second_ap_vector1);
    n2 := crossProduct(ap_vector2, second_ap_vector2);

    common_vector := crossProduct(n1, n2);
    abs_common_vector := vectorModulus(common_vector);
    abs_ap_vector1 := vectorModulus(ap_vector1);
    abs_ap_vector2 := vectorModulus(ap_vector2);

    cosf1 := scalarProduct(ap_vector1, common_vector)/(abs_ap_vector1 * abs_common_vector);
    cosf2 := scalarProduct(ap_vector2, common_vector)/(abs_ap_vector2 * abs_common_vector);

    sinf1 := vectorModulus(crossProduct(ap_vector1, common_vector))/(abs_ap_vector1 * abs_common_vector);

    // sinf1 < 0, если вектор common_vector в 3 или 4 четверти
    if scalarProduct(n1, crossProduct(common_vector, ap_vector1)) > 0 then begin
        sinf1 := -sinf1;
    end;

    sinf2 := vectorModulus(crossProduct(ap_vector2, common_vector))/(abs_ap_vector2 * abs_common_vector);

    // sinf2 < 0, если вектор common_vector в 3 или 4 четверти
    if scalarProduct(n2, crossProduct(common_vector, ap_vector2)) > 0 then begin
        sinf2 := -sinf2;
    end;

    phi := arctan2(sinf1, cosf1);
    p11 := phi; p12 := p11 + Pi;
    r11 := orb1.a * (1 - sqr(orb1.e)) / (1 + orb1.e * cos(p11));
    r12 := orb1.a * (1 - sqr(orb1.e)) / (1 + orb1.e * cos(p12));

    phi := arctan2(sinf2, cosf2);
    p21 := phi; p22 := p21 + Pi;
    r21 := orb2.a * (1 - sqr(orb2.e)) / (1 + orb2.e * cos(p21));
    r22 := orb2.a * (1 - sqr(orb2.e)) / (1 + orb2.e * cos(p22));
end;

// Процедура считывает входные поля с главного окна
procedure readConfigVarsFromScreen(dt, predictionDays : integer; maxDist : Double; form: TForm4);
begin
    dt := StrToInt(form.timeSample.Text);
    maxDist := StrToFloat(form.maxDistance.Text);
    predictionDays := StrToInt(form.predictionTime.Text);
end;

function initSatDatabase(pathToOrbFiles, fileNameORB: string) : array of TSatOrb;
var 
    catDataBase : TCatDataBase;
    isat, d, mon, y, h, m, DAYNR, JD : Integer;
    dmy, hms, maxJD, s : Double;
    satOrbs : array of TSatOrb;
begin
    catDataBase := TCatDataBase.Create(pathToOrbFiles + fileNameORB);
    setLength(satOrbs, catDataBase.catlen);
    maxJD := 0;

    // Считывание орбит фотографа и КО
    for isat := 0 to catDataBase.catlen - 1 do begin
        // Вектор состояния КО
        satOrbs[isat].orbK.a  := catDataBase.cat_a[isat];
        satOrbs[isat].orbK.e  := catDataBase.cat_e[isat];
        satOrbs[isat].orbK.i  := catDataBase.cat_i[isat] * pi / 180;
        satOrbs[isat].orbK.Ap := catDataBase.cat_argp[isat] * pi / 180;
        satOrbs[isat].orbK.RA := catDataBase.cat_dvu[isat] * pi / 180;
        satOrbs[isat].orbK.V  := catDataBase.cat_argsh[isat] * pi / 180 - catDataBase.cat_argp[isat] * pi / 180;
        satOrbs[isat].Kdelta  := catDataBase.cat_bal[isat];

        // Дата, время, юлианская дата, в который измерили вектор состояния КО
        dmy := catDataBase.cat_dmy[isat];
        d := trunc(dmy / 1000000.0);
        mon := trunc((dmy - d * 1000000.0) / 10000);
        Y := round(dmy - trunc(dmy / 10000.0) * 10000);
        JulDay(Y, mon, d, DAYNR, JD);
        hms := catDataBase.cat_hms[isat];
        h := trunc(hms / 10000.0);
        m := trunc((hms - h * 10000.0) / 100);
        s := hms - trunc(hms / 100.0) * 100.0;
        satOrbs[isat].JDDouble := JD + h / 24.0 + m / 1440.0 + s / 86400.0;

        KeplerToGNSK(satOrbs[isat].orbK, satOrbs[isat].orbX);

        // Максимальная юл.дата по базе - для выравнивания всех КО базы к это моменту времени
        if satOrbs[isat].JDDouble > maxJD then begin
            maxJD := satOrbs[isat].JDDouble;
        end;
    end;

    for isat := 0 to length(satOrbs) - 1 do begin
        prognoz_T(satOrbs[isat].orbK.a, satOrbs[isat].orbK.e, satOrbs[isat].orbK.i,
                  satOrbs[isat].orbK.ra, satOrbs[isat].orbK.ap, satOrbs[isat].orbK.v + satOrbs[isat].orbK.ap,
                  0, (maxJD - satOrbs[isat].JDDouble) * 86400,
                  satOrbs[isat].orbX.x, satOrbs[isat].orbX.y, satOrbs[isat].orbX.z,
                  satOrbs[isat].orbX.vx, satOrbs[isat].orbX.vy, satOrbs[isat].orbX.vz);
        GNSKToKepler(satOrbs[isat].orbX, satOrbs[isat].orbK);

        satOrbs[isat].JDDouble := maxJD;
    end;

    result := satOrbs;

    catDataBase.Destroy;
end;

procedure TForm4.Button1Click(Sender: TObject);
var
    isat : integer;
    photoOrbTemp, satOrbTemp: TIOrb;
    closeSats, results, calcAccuracyTest: system.Text;
    orbL : TLagrange;
    start : Double;
    satOrbs : array of TSatOrb; // орбиты КО
    photoOrb: TSatOrb; // орбита фотографа (она же satOrbs[0])
    maxJD, seconds: Double;
    stepNum: integer;
    spaceNum, spaceCount: integer;
    Xstart, Xend: TXYZVxVyVz;
    satMeetingSeries, closeOrbsSeries : array of TPointSeries;
    wasCloseToPhotograph: array of boolean; // запись того, с какими КО было сближение фотографа
    distBetweenSats, maxDist, R: Double;
    dt, satNumber, metSatNum, predictionDays: integer;
    r11, r12, r21, r22 : Double;
    p11, p12, p21, p22 : Double;
begin
    // Создадим файл для записи номеров КО, которые удалось сфотографировать
    system.Assign(closeSats,'..\..\results\Close_satellites_log.txt');
    rewrite(closeSats);
    writeln(closeSats,'Было сближение фотографа со КО: ');

    // Считаем дискрет времени и максимальное расстояние для фотографирования
    readConfigVarsFromScreen(dt, predictionDays, maxDist, Form4);

    // Зададим базу КО
    satOrbs := initSatDatabase('..\..\additional files\', 'orbvnoko150101');

    setLength(wasCloseToPhotograph, length(satOrbs));

    satNumber := length(satOrbs);

    setLength(satMeetingSeries,length(satOrbs));
    setLength(closeOrbsSeries,length(satOrbs));

    // создание трендов, изображающих сближения КО
    for isat := 0 to length(satOrbs) - 1 do begin
        if isat <> 0 then begin
            satMeetingSeries[isat] := TPointSeries.Create(Self);
            satMeetingSeries[isat].ParentChart := Chart1;
            satMeetingSeries[isat].Pointer.Size := 3;
            satMeetingSeries[isat].Pointer.Style := psCircle;
            closeOrbsSeries[isat] := TPointSeries.Create(Self);
            closeOrbsSeries[isat].ParentChart := Chart2;
            closeOrbsSeries[isat].Pointer.Size := 3;
            closeOrbsSeries[isat].Pointer.Style := psCircle;
        end;
    end;

    photoOrb := satOrbs[0];
    stepNum := 1;
    metSatNum := 0;

    // Запуск КО на заданный промежуток времени

    // Цикл по времени
    while (stepNum <= predictionDays * 24 * 3600/dt )  do begin
        writeln(results);
        seconds := stepNum * dt;
        write(results, seconds:7:0, ' ');

        prognoz_T(photoOrb.orbK.a, photoOrb.orbK.e, photoOrb.orbK.i,
                  photoOrb.orbK.ra, photoOrb.orbK.ap, photoOrb.orbK.v + photoOrb.orbK.ap,
                  0, seconds,
                  photoOrb.orbX.x, photoOrb.orbX.y, photoOrb.orbX.z,
                  photoOrb.orbX.vx, photoOrb.orbX.vy, photoOrb.orbX.vz);
        GNSKToKepler(photoOrb.orbX, photoOrbTemp);
        photoOrb.JDDouble := photoOrb.JDDouble + dt/86400;

        // Прогноз движения всех КО на следующие dt секунд
        for isat := 1 to satNumber do begin
            prognoz_T(satOrbs[isat].orbK.a, satOrbs[isat].orbK.e, satOrbs[isat].orbK.i,
                      satOrbs[isat].orbK.ra, satOrbs[isat].orbK.ap, satOrbs[isat].orbK.v + satOrbs[isat].orbK.ap,
                      0, seconds,
                      satOrbs[isat].orbX.x, satOrbs[isat].orbX.y, satOrbs[isat].orbX.z,
                      satOrbs[isat].orbX.vx, satOrbs[isat].orbX.vy, satOrbs[isat].orbX.vz);
            GNSKToKepler(satOrbs[isat].orbX, satOrbTemp);
            satOrbs[isat].JDDouble := satOrbs[isat].JDDouble + dt/86400;

            distBetweenSats := distanceBetweenSatellites(satOrbs[isat].orbX, photoOrb.orbX, dt);

            R := sqrt(sqr(satOrbs[isat].orbX.x - photoOrb.orbX.x) +
                      sqr(satOrbs[isat].orbX.y - photoOrb.orbX.y) +
                      sqr(satOrbs[isat].orbX.z - photoOrb.orbX.z));

            minDistanceBetweenOrbits(photoOrbTemp, satOrbTemp, r11, r12, r21, r22, p11, p12, p21, p22);

            if (abs(r11 - r21) < maxDist) or (abs(r12 - r22) < maxDist) then begin
                closeOrbsSeries[isat].AddXY(seconds, isat, '', clRed);
            end;

            if distBetweenSats <= maxDist then begin
                write(results,'   *');

                satMeetingSeries[isat].AddXY(seconds, isat, '', clGreen);

                if not wasCloseToPhotograph[isat] then begin
                    wasCloseToPhotograph[isat] := True;
                    inc(metSatNum);
                    writeln(closeSats,isat, ' ', satOrbs[isat].JDDouble:7:0);
                end;
            end else write(results,'    ');
        end;

        // Обновляем форму при прогнозе на каждые новые сутки
        if ((stepNum * dt) mod 86400 = 0) then begin
            currentPredictionDay.Caption := 'Текущий день прогноза: ' + IntToStr((stepNum * dt) div 86400)
                 + sLineBreak + 'Обнаружено КО-в: ' + IntToStr(metSatNum);
            Application.ProcessMessages;
//        Self.Refresh;
        end;

        if ((stepNum * dt) div 86400 = predictionDays) then begin
            currentPredictionDay.Caption := currentPredictionDay.Caption + sLineBreak + 'ПРОГНОЗ ОКОНЧЕН';
        end;

        inc(stepNum);
    end;

    system.Close(closeSats);

    while not clearDiagrams.Checked do begin
        Application.ProcessMessages;
    end;

    if clearDiagrams.Checked then begin
        for isat := 0 to length(satOrbs) - 1 do begin
            if isat <> 0 then begin
                satMeetingSeries[isat].Destroy;
                closeOrbsSeries[isat].Destroy;
            end;
        end;
    end;
end;

end.
