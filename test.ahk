#MaxThreadsPerHotkey 2
F12::
    file := FileOpen("coordinates.txt","a")
    toggle:=!toggle
    While toggle{
        MouseGetPos, xpos, ypos
        var := A_TickCount
        file.write("X: " xpos "`tY: " ypos "`t")
        Sleep 16
        file.write("Ticks: " A_TickCount - var "`n")
    }
    file.close()
    
Return