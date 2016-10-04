object Form4: TForm4
  Left = 0
  Top = 0
  Caption = #1044#1080#1072#1075#1088#1072#1084#1084#1099' '#1089#1087#1091#1090#1085#1080#1082#1072'-'#1092#1086#1090#1086#1075#1088#1072#1092#1072
  ClientHeight = 654
  ClientWidth = 1026
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object maxDistanceLabel: TLabel
    Left = 8
    Top = 8
    Width = 118
    Height = 26
    Caption = 'Max '#1088#1072#1089#1089#1090#1086#1103#1085#1080#1077' '#1076#1083#1103' '#1092#1086#1090#1086#1075#1088#1072#1092#1080#1088#1086#1074#1072#1085#1080#1103', '#1082#1084
    Layout = tlCenter
    WordWrap = True
  end
  object satelliteNumberLabel: TLabel
    Left = 152
    Top = 8
    Width = 131
    Height = 26
    Caption = #1063#1080#1089#1083#1086' '#1092#1086#1090#1086#1075#1088#1072#1092#1080#1088#1091#1077#1084#1099#1093' '#1089#1087#1091#1090#1085#1080#1082#1086#1074
    Layout = tlCenter
    WordWrap = True
  end
  object timeSampleLabel: TLabel
    Left = 295
    Top = 8
    Width = 100
    Height = 13
    Caption = #1044#1080#1089#1082#1088#1077#1090' '#1074#1088#1077#1084#1077#1085#1080', '#1089
    Layout = tlCenter
    WordWrap = True
  end
  object predictionTimeLabel: TLabel
    Left = 438
    Top = 8
    Width = 129
    Height = 26
    Caption = #1044#1083#1080#1090#1077#1083#1100#1085#1086#1089#1090#1100' '#1087#1088#1086#1075#1085#1086#1079#1072', '#1076#1085#1077#1081
    Layout = tlCenter
    WordWrap = True
  end
  object currentPredictionDay: TLabel
    Left = 31
    Top = 76
    Width = 135
    Height = 13
    Caption = #1058#1077#1082#1091#1097#1080#1081' '#1076#1077#1085#1100' '#1087#1088#1086#1075#1085#1086#1079#1072': 0'
  end
  object Button1: TButton
    Left = 840
    Top = 33
    Width = 129
    Height = 25
    Caption = #1047#1072#1087#1091#1089#1090#1080#1090#1100' '#1087#1088#1086#1075#1085#1086#1079
    TabOrder = 0
    OnClick = Button1Click
  end
  object Chart1: TChart
    Left = 31
    Top = 136
    Width = 1019
    Height = 242
    Title.Text.Strings = (
      #1044#1080#1072#1075#1088#1072#1084#1084#1072' '#1089#1073#1083#1080#1078#1077#1085#1080#1081' '#1089' '#1092#1086#1090#1086#1075#1088#1072#1092#1086#1084)
    LeftAxis.Automatic = False
    LeftAxis.AutomaticMaximum = False
    LeftAxis.AutomaticMinimum = False
    LeftAxis.Maximum = 51.000000000000000000
    View3D = False
    Zoom.Animated = True
    TabOrder = 1
    DefaultCanvas = 'TGDIPlusCanvas'
    ColorPaletteIndex = 15
    object Series1: TPointSeries
      Title = 'Graph'
      ClickableLine = False
      Pointer.Brush.Gradient.EndColor = 7028779
      Pointer.Gradient.EndColor = 7028779
      Pointer.HorizSize = 1
      Pointer.InflateMargins = True
      Pointer.Style = psCircle
      Pointer.Transparency = 28
      Pointer.VertSize = 1
      Transparency = 28
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
  end
  object maxDistance: TEdit
    Left = 8
    Top = 37
    Width = 138
    Height = 21
    TabOrder = 2
    Text = '50'
  end
  object satelliteNumber: TEdit
    Left = 152
    Top = 37
    Width = 137
    Height = 21
    TabOrder = 3
    Text = '200'
  end
  object timeSample: TEdit
    Left = 295
    Top = 37
    Width = 137
    Height = 21
    TabOrder = 4
    Text = '30'
  end
  object predictionTime: TEdit
    Left = 438
    Top = 37
    Width = 137
    Height = 21
    TabOrder = 5
    Text = '3'
  end
  object Chart2: TChart
    Left = 31
    Top = 384
    Width = 1019
    Height = 249
    Title.Text.Strings = (
      #1044#1080#1072#1075#1088#1072#1084#1084#1072' '#1089#1073#1083#1080#1078#1077#1085#1080#1081' '#1086#1088#1073#1080#1090' '#1087#1086' '#1083#1080#1085#1080#1080' '#1091#1079#1083#1086#1074)
    LeftAxis.Automatic = False
    LeftAxis.AutomaticMaximum = False
    LeftAxis.AutomaticMinimum = False
    LeftAxis.Maximum = 51.000000000000000000
    View3D = False
    Zoom.Animated = True
    TabOrder = 6
    DefaultCanvas = 'TGDIPlusCanvas'
    PrintMargins = (
      15
      41
      15
      41)
    ColorPaletteIndex = 15
    object PointSeries1: TPointSeries
      Title = 'Graph'
      ClickableLine = False
      Pointer.Brush.Gradient.EndColor = 7028779
      Pointer.Gradient.EndColor = 7028779
      Pointer.HorizSize = 1
      Pointer.InflateMargins = True
      Pointer.Style = psCircle
      Pointer.Transparency = 28
      Pointer.VertSize = 1
      Transparency = 28
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {
        00190000000000000000D0814000000000005C82400000000000BC8140000000
        0000F07E400000000000607D400000000000207C400000000000E07F40000000
        0000907F40000000000084824000000000005C82400000000000687F40000000
        000020824000000000002C85400000000000F084400000000000D08140000000
        0000EC834000000000008482400000000000D482400000000000008440000000
        0000BC81400000000000107D400000000000187A4000000000002C8040000000
        0000E481400000000000408540}
    end
  end
  object clearDiagrams: TCheckBox
    Left = 840
    Top = 72
    Width = 146
    Height = 17
    Caption = #1054#1095#1080#1089#1090#1080#1090#1100' '#1076#1080#1072#1075#1088#1072#1084#1084#1099
    TabOrder = 7
  end
end
