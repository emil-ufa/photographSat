program PhotographSat;

uses
  Vcl.Forms,
  Photograph in 'Photograph.pas' {Form4},
  sat_proc in 'sat_proc.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm4, Form4);
  Application.Run;
end.
