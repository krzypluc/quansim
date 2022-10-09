/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Tests for {@link FP64}.
 */
class FP64Test {
    @Test
    void testConsistencyWithDouble() {
        final double v = -5.67e89;
        final Double a = Double.valueOf(v);
        final FP64 b = FP64.of(v);

        Assertions.assertEquals(a.doubleValue(), b.doubleValue());
        Assertions.assertEquals(a.floatValue(), b.floatValue());
        Assertions.assertEquals(a.intValue(), b.intValue());
        Assertions.assertEquals(a.longValue(), b.longValue());
        Assertions.assertEquals(a.byteValue(), b.byteValue());
        Assertions.assertEquals(a.hashCode(), b.hashCode());
    }

    @Test
    void testEquals() {
        final FP64 a = FP64.of(1.23);
        final FP64 b = FP64.of(4.56);

        // Same instance
        Assertions.assertEquals(a, a);
        // Same value
        Assertions.assertEquals(a, FP64.of(a.doubleValue()));
        // Different value
        Assertions.assertNotEquals(a, b);
        // Different object
        Assertions.assertNotEquals(a, new Object());
        Assertions.assertNotEquals(a, null);
    }

    @Test
    void testOne() {
        Assertions.assertEquals(1d, FP64.of(-3.4).one().doubleValue());
    }
    @Test
    void testZero() {
        Assertions.assertEquals(0d, FP64.of(-3.4).zero().doubleValue());
    }

    @Test
    void testSubtract() {
        final double a = 123.4;
        final double b = 5678.9;

        Assertions.assertEquals(a - b, FP64.of(a).subtract(FP64.of(b)).doubleValue());
    }
    @Test
    void testDivide() {
        final double a = 123.4;
        final double b = 5678.9;

        Assertions.assertEquals(a / b, FP64.of(a).divide(FP64.of(b)).doubleValue());
    }

    @Test
    void testMultiplyInt() {
        final double a = 123.4;
        final int n = 3456789;

        Assertions.assertEquals(n * a, FP64.of(a).multiply(n).doubleValue());
    }

    @Test
    void testPowInt() {
        final double a = 123.4;
        final int n = 5;

        Assertions.assertEquals(Math.pow(a, n), FP64.of(a).pow(n).doubleValue());
    }
    @Test
    void testZeroPow() {
        Assertions.assertSame(FP64.of(9876.5).one(), FP64.of(2.3456).pow(0));
    }

    @Test
    void testCompare() {
        Assertions.assertTrue(FP64.of(0).compareTo(FP64.of(-1)) > 0);
        Assertions.assertTrue(FP64.of(1).compareTo(FP64.of(2)) < 0);

        final double v = 123.45;
        Assertions.assertEquals(0, FP64.of(v).compareTo(FP64.of(v)));
    }
}
