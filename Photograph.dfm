object Form4: TForm4
  Left = 0
  Top = 0
  Align = alCustom
  Caption = #1044#1080#1072#1075#1088#1072#1084#1084#1099' '#1089#1087#1091#1090#1085#1080#1082#1072'-'#1092#1086#1090#1086#1075#1088#1072#1092#1072
  ClientHeight = 430
  ClientWidth = 791
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
    Left = 40
    Top = 39
    Width = 118
    Height = 26
    Caption = 'Max '#1088#1072#1089#1089#1090#1086#1103#1085#1080#1077' '#1076#1083#1103' '#1092#1086#1090#1086#1075#1088#1072#1092#1080#1088#1086#1074#1072#1085#1080#1103', '#1082#1084
    Layout = tlCenter
    WordWrap = True
  end
  object timeSampleLabel: TLabel
    Left = 199
    Top = 39
    Width = 100
    Height = 13
    Caption = #1044#1080#1089#1082#1088#1077#1090' '#1074#1088#1077#1084#1077#1085#1080', '#1089
    Layout = tlCenter
    WordWrap = True
  end
  object currentPredictionDay: TLabel
    Left = 63
    Top = 107
    Width = 135
    Height = 13
    Caption = #1058#1077#1082#1091#1097#1080#1081' '#1076#1077#1085#1100' '#1087#1088#1086#1075#1085#1086#1079#1072': 0'
  end
  object spottedSatListLabel: TLabel
    Left = 424
    Top = 168
    Width = 320
    Height = 13
    Caption = #1057#1087#1080#1089#1086#1082' '#1089#1092#1086#1090#1086#1075#1088#1072#1092#1080#1088#1086#1074#1072#1085#1085#1099#1093' '#1050#1054' ('#1076#1077#1085#1100', id '#1050#1054', '#1074#1088#1077#1084#1103' '#1086#1090#1084#1077#1090#1082#1080')'
    WordWrap = True
  end
  object speedStorageLabel: TLabel
    Left = 63
    Top = 148
    Width = 216
    Height = 13
    Caption = #1047#1072#1087#1072#1089' '#1093#1072#1088#1072#1082#1090#1077#1088#1080#1089#1090#1080#1095#1077#1089#1082#1086#1081' '#1089#1082#1086#1088#1086#1089#1090#1080', '#1084'/'#1089': '
  end
  object distLimitLabel: TLabel
    Left = 63
    Top = 215
    Width = 169
    Height = 26
    Caption = 'Max '#1088#1072#1089#1089#1090#1086#1103#1085#1080#1077' '#1084'/'#1091' '#1092#1086#1090#1086#1075#1088#1072#1092#1086#1084' '#1080' '#1050#1054' '#1076#1083#1103' '#1089#1073#1083#1080#1078#1077#1085#1080#1103', '#1082#1084
    WordWrap = True
  end
  object maxSpeedChangeLabel: TLabel
    Left = 241
    Top = 215
    Width = 128
    Height = 26
    Caption = 'Max '#1080#1079#1084#1077#1085#1077#1085#1080#1077' '#1089#1082#1086#1088#1086#1089#1090#1080' '#1087#1088#1080' '#1084#1072#1085#1077#1074#1088#1077', '#1084'/'#1089
    WordWrap = True
  end
  object speedLabel: TLabel
    Left = 241
    Top = 269
    Width = 154
    Height = 15
    Caption = #1047#1072#1087#1072#1089' '#1093#1072#1088'-'#1081' '#1089#1082#1086#1088#1086#1089#1090#1080', '#1084'/'#1089
    WordWrap = True
  end
  object startButton: TButton
    Left = 424
    Top = 34
    Width = 129
    Height = 25
    Caption = #1047#1072#1087#1091#1089#1090#1080#1090#1100' '#1087#1088#1086#1075#1085#1086#1079
    TabOrder = 0
    OnClick = startButtonClick
  end
  object maxDistance: TEdit
    Left = 40
    Top = 68
    Width = 138
    Height = 21
    TabOrder = 1
    Text = '50'
  end
  object timeSample: TEdit
    Left = 199
    Top = 66
    Width = 137
    Height = 21
    TabOrder = 2
    Text = '60'
  end
  object clearDiagrams: TCheckBox
    Left = 424
    Top = 73
    Width = 146
    Height = 17
    Caption = #1054#1095#1080#1089#1090#1080#1090#1100' '#1076#1080#1072#1075#1088#1072#1084#1084#1099
    TabOrder = 3
  end
  object spottedSatList: TListBox
    Left = 424
    Top = 187
    Width = 353
    Height = 201
    ItemHeight = 13
    TabOrder = 4
  end
  object maneuverIsActive: TCheckBox
    Left = 424
    Top = 96
    Width = 146
    Height = 35
    Caption = #1040#1082#1090#1080#1074#1085#1099#1081' '#1088#1077#1078#1080#1084' ('#1089' '#1084#1072#1085#1077#1074#1088#1072#1084#1080')'
    TabOrder = 5
    WordWrap = True
  end
  object distLimit: TEdit
    Left = 63
    Top = 242
    Width = 137
    Height = 21
    TabOrder = 6
    Text = '4000'
  end
  object maxSpeedChange: TEdit
    Left = 241
    Top = 242
    Width = 137
    Height = 21
    TabOrder = 7
    Text = '10'
  end
  object speedStorage: TEdit
    Left = 241
    Top = 290
    Width = 137
    Height = 21
    TabOrder = 8
    Text = '100'
  end
end
