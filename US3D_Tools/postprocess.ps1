$a=$pwd
Get-ChildItem -Path $folderPath | Where { $_.PsIsContainer } | ForEach {
    If ( -Not ($_.Name -match 'Solution')  -and -Not ($_.Name -match 'Grid')){
        echo $_.Name
        cd $_.FullName
      #  Copy-Item "E:\Research\CFD\RAMCII\park90correct2\AddFluxToStagnationLine.m" -Destination ./ -Force
      #  Get-ChildItem -Path $folderPath | Where { $_.PsIsContainer } | ForEach {
      #      Copy-Item "E:\Research\CFD\RAMCII\park90correct2\61km\staganation_line.lay" -Destination $_.FullName -Force
      #  }
       # &matlab -nodisplay -nodesktop -noFigureWindows -nosplash -r LinearizeHeatFlux,quit
        &matlab -nodisplay -nodesktop -noFigureWindows -nosplash -r AddFluxToStagnationLine,quit
    }
}

cd $a