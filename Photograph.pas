﻿unit Photograph;

interface
uses
    Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
    Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, VCLTee.TeEngine,
    VCLTee.Series, Vcl.ExtCtrls, VCLTee.TeeProcs, VCLTee.Chart, Vcl.ComCtrls, System.DateUtils,
    Bal_Types, anomal, astronom, deftorm, numint, atmos, CatDB, precpred, pg_type, math,
    NoradDB, orbint, sgp_h, Lagrange, VCLTee.GanttCh, prognozt, sat_proc, VclTee.TeeGDIPlus;

type
    TForm4 = class(TForm)
    Button1: TButton;
    Chart1: TChart;
    Series1: TPointSeries;
    maxDistance: TEdit;
    maxDistanceLabel: TLabel;
    timeSampleLabel: TLabel;
    timeSample: TEdit;
    predictionTimeLabel: TLabel;
    predictionTime: TEdit;
    clearDiagrams: TCheckBox;
    currentPredictionDay: TLabel;
        {Запускает программу построения диаграммы}
        procedure Button1Click(Sender: TObject);
        private
            { Private declarations }
        public
            { Public declarations }
    end;

var
    Form4: TForm4;

implementation
{$R *.dfm}


procedure TForm4.Button1Click(Sender: TObject);
var
    isat : integer;
    photoOrbTemp, satOrbTemp: TIOrb;
    closeSats, results, calcAccuracyTest: system.Text;
    orbL : TLagrange;
    start : Double;
    satOrbs : TSatOrbBase; // орбиты КО
    photoOrb: TSatOrb; // орбита фотографа (она же satOrbs[0])
    maxJD, seconds: Double;
    stepNum: integer;
    spaceNum, spaceCount: integer;
    Xstart, Xend: TXYZVxVyVz;
    satMeetingSeries : array of TPointSeries;
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
    dt := StrToInt(timeSample.Text);
    maxDist := StrToFloat(maxDistance.Text);
    predictionDays := StrToInt(predictionTime.Text);

    // Зададим базу КО
    initSatDatabase('..\..\additional files\', 'orbvnoko150101', satOrbs);

    setLength(wasCloseToPhotograph, length(satOrbs));

    satNumber := length(satOrbs);

    setLength(satMeetingSeries,length(satOrbs));

    // создание трендов, изображающих сближения КО
    for isat := 0 to length(satOrbs) - 1 do begin
        if isat <> 0 then begin
            satMeetingSeries[isat] := TPointSeries.Create(Self);
            satMeetingSeries[isat].ParentChart := Chart1;
            satMeetingSeries[isat].Pointer.Size := 3;
            satMeetingSeries[isat].Pointer.Style := psCircle;
        end;
    end;

    photoOrb := satOrbs[0];
    stepNum := 1;
    metSatNum := 0;

    // Запуск КО на заданный промежуток времени

    // Цикл по времени
    while (stepNum <= predictionDays * 24 * 3600/dt )  do begin
        seconds := stepNum * dt;

        prognoz_T(photoOrb.orbK.a, photoOrb.orbK.e, photoOrb.orbK.i,
                  photoOrb.orbK.ra, photoOrb.orbK.ap, photoOrb.orbK.v + photoOrb.orbK.ap,
                  0, seconds,
                  photoOrb.orbX.x, photoOrb.orbX.y, photoOrb.orbX.z,
                  photoOrb.orbX.vx, photoOrb.orbX.vy, photoOrb.orbX.vz);
        GNSKToKepler(photoOrb.orbX, photoOrbTemp);
        photoOrb.JDDouble := photoOrb.JDDouble + dt/86400;

        // Прогноз движения всех КО на следующие dt секунд
        for isat := 1 to satNumber - 1 do begin
            if not wasCloseToPhotograph[isat] then begin
                prognoz_T(satOrbs[isat].orbK.a, satOrbs[isat].orbK.e, satOrbs[isat].orbK.i,
                          satOrbs[isat].orbK.ra, satOrbs[isat].orbK.ap, satOrbs[isat].orbK.v + satOrbs[isat].orbK.ap,
                          0, seconds,
                          satOrbs[isat].orbX.x, satOrbs[isat].orbX.y, satOrbs[isat].orbX.z,
                          satOrbs[isat].orbX.vx, satOrbs[isat].orbX.vy, satOrbs[isat].orbX.vz);
                GNSKToKepler(satOrbs[isat].orbX, satOrbTemp);
                satOrbs[isat].JDDouble := satOrbs[isat].JDDouble + dt/86400;

                distBetweenSats := distanceBetweenSatellites(satOrbs[isat].orbX, photoOrb.orbX, dt);

                minDistanceBetweenOrbits(photoOrbTemp, satOrbTemp, r11, r12, r21, r22, p11, p12, p21, p22);

                if distBetweenSats <= maxDist then begin
                    satMeetingSeries[isat].AddXY(seconds, isat, '', clGreen);

                    if not wasCloseToPhotograph[isat] then begin
                        wasCloseToPhotograph[isat] := True;
                        inc(metSatNum);
                    end;
                end;
            end;
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
            end;
        end;
    end;
end;

end.
