#include <SoftwareSerial.h>
SoftwareSerial bluetooth(10, 11);
unsigned long lastSample = 0;

void setup() {
  Serial.begin(115200);
  bluetooth.begin(9600);
  analogReference(DEFAULT);
  pinMode(A0, INPUT);
  lastSample = micros();
}

void loop() {
  unsigned long now = micros();
  if ((unsigned long)(now - lastSample) >= 10000UL) {
    lastSample = now;
    int value = analogRead(A0);
    Serial.println(value);
    bluetooth.print("*Q");
    bluetooth.print(value);
    bluetooth.print('*');
  }
}
// Nominal 100 Hz request interval, not a verified sampling frequency.
// Measure timing and loss before assigning seconds to the phone x-axis.

// Current frame: *Q<integer ADC count>*. TXD->D10; D11->divider->RXD.
// Fixed Q200/Q800 tests passed; this current sketch reads genuine A0 data.
