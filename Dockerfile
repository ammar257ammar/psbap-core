FROM maven:3-jdk-8-alpine as maven

COPY pom.xml .

RUN mvn verify clean --fail-never

COPY src/ ./src/

RUN mvn package 


FROM openjdk:8-jre-alpine

LABEL maintainer "Ammar Ammar <ammar257ammar@gmail.com>"

RUN apk update && apk add bash

COPY --from=maven target/psbap-core-0.0.1-SNAPSHOT-jar-with-dependencies.jar /app/psbap-core.jar

ENTRYPOINT ["java","-jar","/app/psbap-core.jar"]
CMD ["-h"]
