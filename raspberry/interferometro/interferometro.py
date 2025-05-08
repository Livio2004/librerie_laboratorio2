from machine import ADC, Pin
from utime import sleep

ldr = ADC(Pin("GP26"))
led = Pin("GP25", Pin.OUT)
data = []
print("LDR starts reading...")
while True:
    try:
        value = 65535 - ldr.read_u16()  # Read the ADC value
        data.append(value)
        sleep(0.001)  # sleep 0.001 sec

    except KeyboardInterrupt:
        with open('data.txt' , 'w') as f:
            for i in range(len(data)):
                f.write(str(data[i]) + '\n')
        print("finito")
        sleep(0.1)