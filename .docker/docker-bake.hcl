# docker-bake.hcl
variable "ORGANIZATION" {
  default = "aiidateam"
}

variable "REGISTRY" {
  default = "docker.io/"
}

variable "PLATFORMS" {
  default = ["linux/amd64"]
}

variable "QE_VERSION" {
  default = "7.2"
}

function "tags" {
  params = [image]
  result = [
    "${REGISTRY}${ORGANIZATION}/${image}:newly-baked"
  ]
}

group "default" {
  targets = ["aiida-quantumespresso"]
}

target "aiida-quantumespresso" {
  tags = tags("aiida-quantumespresso")
  contexts = {
    src = ".."
  }
  platforms = "${PLATFORMS}"
  args = {
    "QE_VERSION" = "${QE_VERSION}"
  }
}