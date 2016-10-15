unit sat_proc;

interface

uses Lagrange, bal_types, math, CatDB, astronom, System.SysUtils, prognozt;

type
   TSatOrb = record
        id : Integer;
        JDDouble : Double;
        Kdelta : Double;
        orbX : TXYZVxVyVz;
        case integer of
            0 : (orbK: TIOrb;);
            1 : (items : array [0..5] of Double;);
        end;

    TSatMeeting = record
        sat : TSatOrb; // ��, � ������� ����� ����������� ��� ������ �������
        T : Double; // ������ ������� (�������), � ������� ����� ��������� �������
        dv : Double; // ����� ���������� �������� ����� ������� ��������� ��� �������
    end;
    TSatOrbBase = array of TSatOrb;

    function distanceBetweenSatellites(x1, x2 : TXYZVxVyVz; dt : Double): Double;
    procedure minDistanceBetweenOrbits(orb1, orb2 : TIOrb; var r11, r12, r21, r22, p11, p12, p21, p22 : Double);
    procedure initSatDatabase(pathToOrbFiles, fileNameORB : string; var satOrbs : TSatOrbBase);
    function dateForLog() : String;
    function calculateEa(orbK : TIOrb) : double;
    function timeBetweenDifferentV(orbK : TIorb; v1, v2 : double) : double;
    function setKeplerOrb(a, e, i, ap, ra, v : double) : TIorb;
    procedure makeManeuver(dV : double; var orbK : TIorb);

    // ������� ��� �������
    function getCorrectAngle(phi : double) : double;
    function getDistDiff(dv : double; dt, T: integer; photoOrb, satOrb : TSatOrb) : double;
    function getManeuverSpeedChange(dvMax, maxDistDiff, dt, T : Integer; photoOrb, satOrb : TSatOrb) : double;
implementation

// ������� ���������� ����� ���������� � �� � �������� ������ ������� �
// (� ��������) ��� ���������� ���������� ����������� ������� �����������
// �������� dv (�/�). ����� dt - ������� ������� ��� �������� �������� ���������.
//
function getDistDiff(dv : double; dt, T: integer; photoOrb, satOrb : TSatOrb) : double;
var
    orbK : TIOrb;
    p, s : TXYZVxVyVz;
begin
    orbK := photoOrb.orbK;
    makeManeuver(dv, orbK);

    prognoz_T(orbK.a, orbK.e, orbK.i,
              orbK.ra, orbK.ap, orbK.v + orbK.ap,
              0, T,
              p.x, p.y, p.z, p.vx, p.vy, p.vz);

    orbK := satOrb.orbK;
    prognoz_T(orbK.a, orbK.e, orbK.i,
              orbK.ra, orbK.ap, orbK.v + orbK.ap,
              0, T,
              s.x, s.y, s.z, s.vx, s.vy, s.vz);
    result := distanceBetweenSatellites(p, s, dt);
end;

function getCorrectAngle(phi : double) : double;
begin
    if phi < -pi then begin
        while (phi < -pi) do begin
            phi := phi + 2 * pi;
        end;
    end;

    if phi >= pi then begin
        while (phi >= pi) do begin
            phi := phi - 2 * pi;
        end;
    end;

    result := phi;
end;

// ������� ��������� ����������� ���������� �� �������� ��� ���������,
// ����� ���������� � �� ����� ����� T ������ �� ���������� ������� (0 ������).
// ��������������, ��� � ������ T �� ��������� �� ����� ����� (��� ����� ������ � ���)
//
// ���������:
//      dt - ������� ������� ��� �������� �������� ��������� (� ��������)
//      T - ������ ������� (� ��������), �� ������� ������������ ��������� ��������� � ��
//      photoOrb, satOrb - ������� ��������� ��������� � �� ��������������
//
function getManeuverSpeedChange(dvMax, maxDistDiff, dt, T : Integer; photoOrb, satOrb : TSatOrb) : double;
var
    dv, distDiff : double;
begin
    dv := -dvMax;
    while dv <= dvMax do begin
        distDiff := getDistDiff(dv, dt, T, photoOrb, satOrb);
        if distDiff < maxDistDiff then break;

        dv := dv + 0.01;
    end;

    if (dv = dvMax) and (distDiff > maxDistDiff) then result := 0
    else result := dv;
end;

// ��������� ��������� ���������� ������ ��
//
// ������� ���������:
// dV - �������� ���������� �� �������� � �/� (����� ���� <0)
// orbK - ������ ��������� �� � ���������� �����������
procedure makeManeuver(dV: double; var orbK: TIorb);
var
    orbX: TXYZVxVyVz;
    v: double; // �������� ��
begin
    KeplerToGNSK(orbK, orbX);
    v := vectorModulus(orbX.vel);
    orbX.vx := orbX.vx * (1 + dV * 1e-3/v);
    orbX.vy := orbX.vy * (1 + dV * 1e-3/v);
    orbX.vz := orbX.vz * (1 + dV * 1e-3/v);
    GNSKToKepler(orbX, orbK);
end;

// ������� ������������ ���������� ����� ����� ��, ��������� ������� ����������
// ������ dt ������. ���� �� ����� dt ����������� ����������� ���������� �����
// �� (���������� �� ������ �� ��� �� ������, �� ������� ����� ������),
// ����� ������� ���������� ��� ����������� ����������.
//
// ������� ���������:
//      x1, x2 - ������� ��������� ������� ��
//      dt     - ��������� �������� ��������� ��������� ��
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

// ������� ������������ ����������� ���������� ����� ��������
// ���� �� �� ����� ����� (��) - ������ ����������� ���������� �����.
//
// ������� ���������:
//      orb1, orb2 - ������� ��������� ���� �� (� ��������� �������)
//
// �������� ���������:
//      r11, r12, r21, r22 - ���������� �� �������������� ������ �� ����� ����������� ����� � ��
//      p11, p12, p21, p22 - ��������� ����� ����������� ����� � ��
procedure minDistanceBetweenOrbits(orb1, orb2 : TIOrb; var r11, r12, r21, r22, p11, p12, p21, p22: Double);
var
    X : TXYZVxVyVz; // ������������� ������ � ���������� �����������
    orbK: TIOrb; // ������������� �������� �������
    ap_vector1, ap_vector2 : TXYZ; // ������-������� �������� �����
    common_vector : TXYZ; // ������������ ������ ����� �����
    second_ap_vector1, second_ap_vector2: TXYZ; // �������, �������� �������������� ��������� 1� � 2� ������
    n1, n2: TXYZ; // ������� � ���������� �����
    cosf1, sinf1, cosf2, sinf2 : Double; // cos � sin �����, �������� ����� ����������� ����� � ������ �����
    abs_common_vector, abs_ap_vector1, abs_ap_vector2 : Double; // ������ ������������ ��������
    phi : Double; // ������������� ����
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

    // sinf1 < 0, ���� ������ common_vector � 3 ��� 4 ��������
    if scalarProduct(n1, crossProduct(common_vector, ap_vector1)) > 0 then begin
        sinf1 := -sinf1;
    end;

    sinf2 := vectorModulus(crossProduct(ap_vector2, common_vector))/(abs_ap_vector2 * abs_common_vector);

    // sinf2 < 0, ���� ������ common_vector � 3 ��� 4 ��������
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

procedure initSatDatabase(pathToOrbFiles, fileNameORB: string; var satOrbs : TSatOrbBase);
var
    catDataBase : TCatDataBase;
    isat, d, mon, y, h, m, DAYNR, JD : Integer;
    dmy, hms, maxJD, s : Double;
    catLen : Integer;
    satOrbsTemp : TSatOrbBase;
begin
    catDataBase := TCatDataBase.Create(pathToOrbFiles + fileNameORB);
    catLen := catDataBase.catlen;
    setLength(satOrbsTemp, catLen);
    setLength(satOrbs, catLen);
    maxJD := 0;

    // ���������� ����� ��������� � ��
    for isat := 0 to catDataBase.catlen - 1 do begin
        // ������ ��������� ��
        satOrbsTemp[isat].orbK.a  := catDataBase.cat_a[isat];
        satOrbsTemp[isat].orbK.e  := catDataBase.cat_e[isat];
        satOrbsTemp[isat].orbK.i  := catDataBase.cat_i[isat] * pi / 180;
        satOrbsTemp[isat].orbK.Ap := catDataBase.cat_argp[isat] * pi / 180;
        satOrbsTemp[isat].orbK.RA := catDataBase.cat_dvu[isat] * pi / 180;
        satOrbsTemp[isat].orbK.V  := catDataBase.cat_argsh[isat] * pi / 180 - catDataBase.cat_argp[isat] * pi / 180;
        satOrbsTemp[isat].Kdelta  := catDataBase.cat_bal[isat];
        satOrbsTemp[isat].id := isat;

        // ����, �����, ��������� ����, � ������� �������� ������ ��������� ��
        dmy := catDataBase.cat_dmy[isat];
        d := trunc(dmy / 1000000.0);
        mon := trunc((dmy - d * 1000000.0) / 10000);
        Y := round(dmy - trunc(dmy / 10000.0) * 10000);
        JulDay(Y, mon, d, DAYNR, JD);
        hms := catDataBase.cat_hms[isat];
        h := trunc(hms / 10000.0);
        m := trunc((hms - h * 10000.0) / 100);
        s := hms - trunc(hms / 100.0) * 100.0;
        satOrbsTemp[isat].JDDouble := JD + h / 24.0 + m / 1440.0 + s / 86400.0;

        KeplerToGNSK(satOrbsTemp[isat].orbK, satOrbsTemp[isat].orbX);

        // ������������ ��.���� �� ���� - ��� ������������ ���� �� ���� � ��� ������� �������
        if satOrbsTemp[isat].JDDouble > maxJD then begin
            maxJD := satOrbsTemp[isat].JDDouble;
        end;

        satOrbs := satOrbsTemp;
    end;

    for isat := 0 to length(satOrbsTemp) - 1 do begin
        prognoz_T(satOrbsTemp[isat].orbK.a, satOrbsTemp[isat].orbK.e, satOrbsTemp[isat].orbK.i,
                  satOrbsTemp[isat].orbK.ra, satOrbsTemp[isat].orbK.ap, satOrbsTemp[isat].orbK.v + satOrbsTemp[isat].orbK.ap,
                  0, (maxJD - satOrbsTemp[isat].JDDouble) * 86400,
                  satOrbsTemp[isat].orbX.x, satOrbsTemp[isat].orbX.y, satOrbsTemp[isat].orbX.z,
                  satOrbsTemp[isat].orbX.vx, satOrbsTemp[isat].orbX.vy, satOrbsTemp[isat].orbX.vz);
        GNSKToKepler(satOrbsTemp[isat].orbX, satOrbsTemp[isat].orbK);

        satOrbsTemp[isat].JDDouble := maxJD;
    end;

    catDataBase.Destroy;
end;

function dateForLog() : String;
var
    dttm : TDateTime;
    year, month, day, hour, min, sec, msec : word;
begin
    dttm := Now;
    DecodeDate(dttm, Year, Month, Day);
    DecodeTime(dttm, Hour, Min, Sec, MSec);


    result := IntToStr(year) + IntToStr(month) + IntToStr(day)
              + '_' + IntToStr(hour) + IntToStr(min) + IntToStr(sec);
end;

// ������� ������������ ��������������� �������� �� �������� �������� �
// ���������������
function calculateEa(orbK: TIOrb) : double;
var
    e, v, cosEa, sinEa: double;
begin
    v := orbK.v;
    e := orbK.e;
    cosEa := (cos(v) + e)/(1 + e * cos(v));
    sinEa := sin(v) * sqrt(1 - e*e)/(1 + e * cos(v));
    result := ArcTan2(sinEa, cosEa);
end;

// ������� ������������ ����� ������ �������� �� v = v1 �� v = v2
function timeBetweenDifferentV(orbK: TIorb; v1, v2: double) : double;
var
    T : double;
    Ea1, Ea2, e: double;
begin
    while (v1 > 2 * pi) do v1 := v1 - 2*pi;
    while (v1 < -2 * pi) do v1 := v1 + 2*pi;
    while (v2 > 2 * pi) do v2 := v2 - 2*pi;
    while (v2 < -2 * pi) do v2 := v2 + 2*pi;

    T := 2 * pi * sqrt(orbK.a*orbK.a*orbK.a/GM_km);
    e := orbK.e;

    orbK.v := v1;
    Ea1 := calculateEa(orbK);
    orbK.v := v2;
    Ea2 := calculateEa(orbK);

    if v2 >= v1 then begin
        result := T/(2*pi)*((Ea2 - e*sin(Ea2)) - (Ea1 - e*sin(Ea1)));
    end else begin
        result := T - abs(T/(2*pi)*((Ea2 - e*sin(Ea2)) - (Ea1 - e*sin(Ea1))));
    end;
end;

// ������� ��� ������� ���������� ������
function setKeplerOrb(a, e, i, ap, ra, v : double) : TIorb;
begin
    result.a := a;
    result.e := e;
    result.i := i;
    result.ap := ap;
    result.ra := ra;
    result.v := v;
end;

end.
