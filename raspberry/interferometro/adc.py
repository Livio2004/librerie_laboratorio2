from machine import ADC, Pin
from utime import sleep

pin = Pin("GP22", Pin.OUT)
adc = ADC(pin('GP26'))

print("ADC starts reading...")
while True:
    try:
        value = adc.read_u16()  # Read the ADC value
        print("ADC Value:", value)
        sleep(1)  # sleep 1sec
    except KeyboardInterrupt:
        break
pin.off()
print("Finished.")