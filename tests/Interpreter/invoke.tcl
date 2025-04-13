

model basic 1 1
uniaxialMaterial Elastic 1 1 2


invoke UniaxialMaterial 1 {

    foreach strain [linspace 0.0 0.01 20] {

      strain $strain -commit

      puts "$strain [stress]"

    }
}


