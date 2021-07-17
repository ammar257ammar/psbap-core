FROM maven:3-jdk-8-alpine as maven

COPY pom.xml .

RUN mvn verify clean --fail-never

COPY src/ ./src/

RUN mvn package 


FROM ubuntu:bionic-20191029

LABEL maintainer "Ammar Ammar <ammar257ammar@gmail.com>"

RUN apt-get update && \
	apt-get install -y wget openjdk-8-jdk && \
	apt-get install -y ant && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/cache/oracle-jdk8-installer

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

COPY --from=maven target/psnpbind-core-1.0-Stable-jar-with-dependencies.jar /app/psnpbind-core.jar

ENTRYPOINT ["java","-jar","/app/psnpbind-core.jar"]
CMD ["-h"]
